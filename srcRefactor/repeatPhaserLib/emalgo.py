'''
sudo pip install python-Levenshtein

Iteratively compute the EM estimate of the solution

Input : RList_Double.fasta, contigLeft.json, contigRight.json, intermediate.fasta
Output : score, matching, betterInteriorToFlank.fasta 

High level requirement : 
1) Produce the matching of contigs 
2) Produce interior of reads
3) Provide a level of confidence

'''
from srcRefactor.repeatPhaserLib.finisherSCCoreLib import IORobot
from srcRefactor.repeatPhaserLib.finisherSCCoreLib import alignerRobot
from srcRefactor.repeatPhaserLib.finisherSCCoreLib import graphLib
from srcRefactor.repeatPhaserLib.finisherSCCoreLib import houseKeeper
from srcRefactor.repeatPhaserLib import abunGraphLib
from srcRefactor.repeatPhaserLib import abunHouseKeeper
from operator import itemgetter
from itertools import groupby
from itertools import permutations

from Bio.Seq import Seq
import Bio.Align.Applications
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import os 
from Bio.Align import AlignInfo
from Bio import pairwise2
import json
import Levenshtein
import random 
import math
import time
import numpy as np
from multiprocessing import Pool


from srcRefactor.dataGenLib import dataGenLib

def readInJSON(folderName, filename):
	json_data = open(folderName + filename, 'r')
	dataItem = json.load(json_data)
	return dataItem

def preparation(folderName):
	'''
	Prepare RList.fasta, contigLeft.json, contigRight.json, intermediate.fasta
	This step will not be needed in production as it should be automatically given or will follow a different logic of generation	
	'''

	CLeftList , CRightList = [] , [] 
	RList = [] 
	templateList = [] 

	contigReadGraph = "phaseStringGraph1"
	G = graphLib.seqGraph(0)
	G.loadFromFile(folderName, contigReadGraph)

	lenDic = IORobot.obtainLength(folderName, "improved3_Double.fasta")
	N1 = len(lenDic)

	kthres, edgeThres = 3, 1

	G = graphLib.seqGraph(0)
	G.loadFromFile(folderName, contigReadGraph)

	if True:
		adj = [[] for i in range(N1)]

		for i in range(N1): 
			tmpList = abunGraphLib.findAllReachable(i, N1, G)

			for j in tmpList:
				if len(abunGraphLib.findAllPathK(i,j,G,kthres)) >= edgeThres:
					adj[i].append(j) 
		
		print adj

	if True:
		CLeftList , CRightList = [0, 6] , [4, 8] 
		RList = [] 
		templateList = []

		numberOfFiles = houseKeeper.globalParallelFileNum
		dataList = []
		for i in range(1, 1+numberOfFiles): 
			if i < 10:
				indexOfMum = "0" + str(i)
			else:
				indexOfMum = str(i)

			dataList = dataList+ alignerRobot.extractMumData(folderName, "outAbun"+ str(indexOfMum)+"Out")


		middleList =  [2]

		CLeftNameList, CRightNameList, middleNameList \
			= [abunHouseKeeper.parseIDToName(i, 'C', 0) for i in CLeftList]  \
			, [abunHouseKeeper.parseIDToName(i, 'C', 0) for i in CRightList] \
			, [abunHouseKeeper.parseIDToName(i, 'C', 0) for i in middleList]
		
		dataList.sort(key = itemgetter(-2))
		for key, items in groupby(dataList, itemgetter(-2)): 
			print key
			if int(key[5:]) == 1 :
				for eachitem in items:
					RList.append(eachitem[-1])
					#print eachitem[-4]

		#print "len(RList) : ", len(RList)
		RList = abunHouseKeeper.getDistinct(RList)
		print "len(RList) : ", len(RList)
		lenDic = IORobot.obtainLength(folderName, "improved3.fasta")
		print lenDic["Segkk1"]
		# print RList
		IORobot.putListToFileO(folderName, "raw_reads.fasta", "RList", RList)

		ctgList= ["Contig0_p", "Contig3_p"]
		with open(folderName +  "contigLeft.json", 'w') as outfile:
			json.dump(ctgList, outfile)
		

		ctgList= ["Contig2_p", "Contig4_p"]
		with open(folderName +  "contigRight.json", 'w') as outfile:
			json.dump(ctgList, outfile)
		
		contigDic = IORobot.loadContigsFromFile(folderName, "improved3_Double.fasta")

		#addNoise(contigDic["Contig1_p"])
		### no noise 
		# IORobot.writeSegOut([contigDic["Contig1_p"]], folderName, "intermediate.fasta")
		### with noise
		noisyIntermediate = dataGenLib.createANoisyRead(len(contigDic["Contig1_p"]), 0.01, contigDic["Contig1_p"])
		IORobot.writeSegOut([noisyIntermediate], folderName, "intermediate.fasta")
		IORobot.writeSegOut([contigDic["Contig1_p"]], folderName, "intermediateNoiseless.fasta")

def computeReadAssociation(folderName, prevIteration, constants, isDebug, mummerLink):
	'''
	Input :  prevIteration = [lambdas, templates, score] , constants = [basesMappedToEachContig]
	Output : readMatching = [read2templateDic, template2readDic]
	Algorithm : 
	1) Align reads to the templates
	2) Assign them to templates. Tie break by lambda. 
	3) Format return as Dictionary and return 
	'''
	IORobot.writeSegOut(prevIteration[1], folderName, "templates.fasta")
	#assert(False)
	readMatching = findAnchors(folderName, prevIteration, isDebug, mummerLink)
	return readMatching

def computeLambdaList(folderName, prevIteration, readMatching, constants):
	'''
	Input : prevIteration = [lambdas, templates, score] , constants = [basesMappedToEachContig], readMatching = [read2templateDic, template2readDic]
	Output : lambdas 
	Algorithm : 
	1) Use read2templateDic and basesMappedToEachContig to compute the new lambda
	2) Format and return 

	### ideally, length is not used, but only the number... 
	'''
	lambdas = []
	lenDic = IORobot.obtainLength(folderName , "RList_Double.fasta")
	template2readDic = readMatching[1] 

	tmp = [0 for i in range(len(constants[0]))]
	for i in range(len(constants[0])):
		tmp[i] = constants[0][i] 
		for eachitem in template2readDic["Segkk"+ str(i)]:
			tmp[i] += lenDic[eachitem[-2]]

	total = np.sum(tmp)
	for i in range(len(tmp)):
		lambdas.append(tmp[i]*1.0/total)


	# print prevIteration[0]
	# print lambdas
	
	return lambdas 

def computeInteriorList(folderName , prevIteration, readMatching, constants, isDebug, mummerLink):
	'''
	Input : prevIteration = [lambdas, templates, score] , constants = [basesMappedToEachContig], readMatching = [read2templateDic, template2readDic]
	Output : interiors = [newTemplates]
	Algorithm : 
	1) Use read2templateDic, templates, and RList.fasta 
		a)  Set up bins on templates
		b)  Align reads segments to bins 
		c)	Consensus within bins
		d)	Concat the results by back processing the edits
	2) Format and return the interiors
	'''
	if not  isDebug:
		interiors = chopUpReads(folderName, mummerLink)
	else:
		interiors  = IORobot.readContigsFromFile(folderName, "interiors.fasta")
		

	return interiors 

def concatAndFormat(folderName):
	assert(True)

def findAnchors(folderName, prevIteration, isDebug, mummerLink):
	'''
	Input: IORobot.writeSegOut(ctgList, folderName, "templates.fasta")
		   IORobot.putListToFileO(folderName, "raw_reads.fasta", "RList", RList)
    Output : the assignmentDic and lookUpDic

	'''

	if not  isDebug:			
		alignerRobot.useMummerAlign(mummerLink, folderName, "templateAnchor", "RList_Double.fasta", "templates.fasta", False, "", False)
		dataList = alignerRobot.extractMumData(folderName, "templateAnchor" +"Out")
		
		lenDicReads = IORobot.obtainLength(folderName, "RList_Double.fasta")
		
		lenDicTemplates = IORobot.obtainLength(folderName, "templates.fasta")

		templatesDic = IORobot.loadContigsFromFile(folderName, "templates.fasta")
		readsDic = IORobot.loadContigsFromFile(folderName, "RList_Double.fasta")


		#print templatesDic["Segkk0"][1144:1144+ 50]
		#print readsDic["Segkk11098"][47:47+50]
		#print Levenshtein.distance(templatesDic["Segkk0"][1144:1144+ 50], readsDic["Segkk11098"][47:47+50])
		#assert(False)

		read2templateDic = {}
		template2readDic = {}

		for eachitem in lenDicTemplates:
			template2readDic[eachitem] = []

		for eachitem in lenDicReads	:
			read2templateDic[eachitem] = []
		
		thres = 30
		dataList.sort(key = itemgetter(-2))
		for key, items in groupby(dataList, itemgetter(-2)):
			L = lenDicReads[key]
			tmpList = []
			for eachitem in items:
				if eachitem[4] > L - thres and eachitem[2] < eachitem[3]:
					tmpList.append(eachitem)

			if len(tmpList) >= 1: 
				returnItem = resolveCompetingTemplates(folderName, tmpList, key, templatesDic, readsDic, prevIteration[0])
				readName, templateName = returnItem[-2], returnItem[-1]
				read2templateDic[readName].append(returnItem[-1])
				template2readDic[templateName].append(returnItem)


		with open(folderName +  "read2templateDic.json", 'w') as outfile:
			json.dump(read2templateDic, outfile)

		with open(folderName +  "template2readDic.json", 'w') as outfile:
			json.dump(template2readDic, outfile)

		print len(dataList), len(lenDicReads), len(lenDicTemplates), len(read2templateDic)
		# assert(False)


	else:
		read2templateDic = readInJSON(folderName, "read2templateDic.json")
		template2readDic = readInJSON(folderName, "template2readDic.json")


	return [read2templateDic, template2readDic]

def resolveCompetingTemplates(folderName, competingTemplateList, readName,  templatesDic, readsDic, lambdas):
	returnItem = []

	tmpScore = -1

	if len(competingTemplateList) > 1:
		ct = 0
		for eachitem in competingTemplateList:
			tmpTemplateName = eachitem[-1]
			start, end = eachitem[2], eachitem[3]
			score = Levenshtein.distance(readsDic[readName],templatesDic[tmpTemplateName][start:end])

			if score > tmpScore:
				tmpScore = score
				returnItem = eachitem
				ct = lambdas[int(tmpTemplateName[5:])]
			elif score == tmpScore:
				coinFlip = random.random()
				ct += lambdas[int(tmpTemplateName[5:])]
				if coinFlip < lambdas[int(tmpTemplateName[5:])]*1.0/ct: 
					returnItem = eachitem

		return returnItem

	elif len(competingTemplateList) == 1:
		return competingTemplateList[0]
	
	else:
		return returnItem

class consensusBins(object):
	def __init__(self, index, begin, end):
		self.index = index
		self.begin = begin
		self.end = end
		self.readList = [] 

	def addToReadList(self, item):
		self.readList.append(item)
	def printInfo(self):
		print self.index, self.begin, self.end , self.readList
	def fetchReadInfo(self, readName, prevEnd):
		thres = 10

		for eachitem in self.readList:
			if eachitem[0] == readName and abs(eachitem[1] - prevEnd) < thres :
				return eachitem

		return None

def chopUpReads(folderName, mummerLink):
	print "chopUpReads"
	interiors = []
	
	### Initializtion 
	read2templateDic = readInJSON(folderName, "read2templateDic.json")
	template2readDic = readInJSON(folderName, "template2readDic.json")

	for eachitem in template2readDic:
		print "Length", eachitem, len(template2readDic[eachitem])
	#assert(False)
		
	dataList = alignerRobot.extractMumData(folderName, "templateAnchor" +"Out")
	lenDicTemplates = IORobot.obtainLength(folderName, "templates.fasta")
	templatesDic = IORobot.loadContigsFromFile(folderName, "templates.fasta")
	readsDic = IORobot.loadContigsFromFile(folderName, "RList_Double.fasta")
	dataList.sort(key = itemgetter(-1))



	### Set up bins	
	templateBeginEndDic = {}
	ell = 50

	for key, items in groupby(dataList, itemgetter(-1)):
		begin, end = 10**9 ,-1  
		for eachitem in items:
			if eachitem[2] <  begin:
				begin = eachitem[2]

			if eachitem[3] > end:
				end = eachitem[3]

		templateBeginEndDic[key] = [begin, end]

	print templateBeginEndDic
	
	GTDic = IORobot.loadContigsFromFile(folderName, "GTDic.fasta")
		
	for i in range(len(lenDicTemplates)):
		nameOfTemplate = "Segkk"+str(i)
		begin, end = templateBeginEndDic[nameOfTemplate] 
		numberOfBins = int(math.ceil((end - begin)*1.0/ell))
		print "numberOfBins", numberOfBins
		bins = [consensusBins(j, begin + j*ell, min(begin +ell*(j+1), end)) for j in range(numberOfBins)]

		temp2readAlignDic = loadAlignment(folderName, nameOfTemplate, template2readDic[nameOfTemplate],mummerLink)
		#for eachdebug in temp2readAlignDic:
		#	print temp2readAlignDic[eachdebug]
		
		### Align reads to bins
		for eachalign in template2readDic[nameOfTemplate]:
			templateStart, templateEnd, readStart, readEnd = eachalign[2], eachalign[3], eachalign[0], eachalign[1]
			readName = eachalign[-2]
			#print readName, templateStart, templateEnd, readStart, readEnd , nameOfTemplate
			#assert(False)

			indexOfBin = min( math.ceil((templateStart - begin )*1.0/ ell), numberOfBins-1)
			indexOfBin =  int(indexOfBin)

			while bins[indexOfBin].end < templateEnd:
				binStart, binEnd = bins[indexOfBin].begin, bins[indexOfBin].end
				readSegStart, readSegEnd = temp2readAlignDic[readName][binStart], temp2readAlignDic[readName][binEnd] 
				bins[indexOfBin].addToReadList([readName, readSegStart, readSegEnd])
				indexOfBin += 1	 

		timestart = time.time() 
		returnString = localConsensus(folderName, bins, readsDic, templatesDic, nameOfTemplate)
		print time.time() - timestart
		#assert(False)
		interiors.append(returnString)


		print "TemplateDist : ",  Levenshtein.distance(templatesDic[nameOfTemplate], GTDic[nameOfTemplate])
		print "CleanedDist : ", Levenshtein.distance(returnString, GTDic[nameOfTemplate])
		#print returnString[10257-3:10257+3], GTDic[nameOfTemplate][10259-3:10259+3]
		#for eachedit in Levenshtein.editops(returnString, GTDic[nameOfTemplate]):
		#	print eachedit
		# assert(False)
		
	IORobot.writeSegOut(interiors, folderName, "interiors.fasta")	
	return interiors

def loadAlignment(folderName, nameOfTemplate, associatedReadAlignList,mummerLink):


	temp2readAlignDic = {}
	for eachitem in associatedReadAlignList:
		readName = 	eachitem[-2]
		dic = parseAlignment(folderName, nameOfTemplate, readName, mummerLink)
		temp2readAlignDic[readName] = dic

	#assert(False)

	return temp2readAlignDic

def parseAlignment(folderName, nameOfTemplate, readName, mummerLink):
	# print "show-aligns [options] <delta file> <IdR> <IdQ>"
	# /Users/kakitlam/Desktop/experimentBench/MUMmer3.23/show-aligns -w 10000  dataFolder/templateAnchor.delta Segkk48988 Segkk1
	# [1, 5986, 6734, 12733, 5986, 6000, 97.72, u'Segkk26449', u'Segkk0']

	mummerPath = mummerLink
	command = mummerPath + "show-aligns -w 100000 " + folderName + "templateAnchor.delta " + readName + " " + nameOfTemplate + "  | head -10 | tail -2   > " + folderName + "alignmentFile"
	os.system( command)

	f = open(folderName + "alignmentFile", 'r')
	tmp = f.readline().rstrip()
	dataList0 = tmp.split()
	tmp = f.readline()
	dataList1 = tmp.split()
	f.close()
	
	readStart , seqRead= int(dataList0[0]), dataList0[1]
	templateStart, seqTemplate = int(dataList1[0]), dataList1[1]
	
	dic = {}
	readCtr , templateCtr = 0 ,0 
	for i in range(len(seqRead)):
		dic[templateStart + templateCtr] = readStart + readCtr 
		### logic for the offset here...
		if seqRead[i] != '.':
			readCtr += 1

		if seqTemplate[i] != '.':
			templateCtr += 1

	return dic

def calculate(func, args):
	func(*args)

def calculatestar(args):
	calculate(*args)

def localConsensus(folderName, bins, readsDic,templatesDic, nameOfTemplate):
	print "localConsensus"

	totalChangeList = []
	if False:
		for k in range(len(bins)-1):
			if len(bins[k].readList) > 1:
				segList = []		
				segList.append(templatesDic[nameOfTemplate][bins[k].begin - 1 :bins[k+1].end - 1 ])

				for eachseg in bins[k].readList:
					readName, start, prevEnd  = eachseg[0], eachseg[1], eachseg[2]
					readInfo = bins[k+1].fetchReadInfo(readName, prevEnd) 
					if readInfo != None:
						end = readInfo[-1]
						segList.append(readsDic[readName][start- 1 :end- 1 ])

				if len(segList) > 2:
					modiList=  clustalw2MSA(folderName, segList, bins[k].begin)
				else:
					modiList = []

				totalChangeList = totalChangeList + modiList

	else:
		print "Parallel"
		workerList = []
		nProc = houseKeeper.globalParallel

		for k in range(len(bins)-1):
			if len(bins[k].readList) > 1:
				segList = []		
				segList.append(templatesDic[nameOfTemplate][bins[k].begin - 1 :bins[k+1].end - 1 ])

				for eachseg in bins[k].readList:
					readName, start, prevEnd  = eachseg[0], eachseg[1], eachseg[2]
					readInfo = bins[k+1].fetchReadInfo(readName, prevEnd) 
					if readInfo != None:
						end = readInfo[-1]
						segList.append(readsDic[readName][start- 1 :end- 1 ])

				if len(segList) > 2:
					workerList.append([folderName, segList, bins[k].begin])

		p = Pool(processes=nProc)
		results = []

		for eachitem in workerList:
			folderName, segList, startingIndex= eachitem
			results.append(p.apply_async(clustalw2MSA, args=(folderName, segList, startingIndex)))

		for itemkk in results:
			totalChangeList += itemkk.get() 
		
		p.close()

	totalChangeList.sort(reverse = True)
	sequence = templatesDic[nameOfTemplate]
	characterList = ['S']

	for eachitem in sequence:
		characterList.append(eachitem)

	for key , items in groupby(totalChangeList):
		#print key
		loc, typeOfEdit = key[0] , key[1]
		if typeOfEdit[0] == 'd':
			characterList.pop(loc)
		elif typeOfEdit[0] == 's':
			infoArr = typeOfEdit.split('_')
			characterList[loc] = infoArr[-1]
		elif typeOfEdit[0] == 'i':
			infoArr = typeOfEdit.split('_')
			characterList.insert(loc, infoArr[-1])

	characterList.pop(0)
	#assert(False)
	returnStr = ""
	for i in characterList:
		returnStr += i

	return returnStr

def clustalw2MSA(folderName, segList, startingIndex):
	filename =  "seq" + str(startingIndex)
	lenMax = min(11, len(segList))

	'''
	indexNum = 10259
	if  False and  startingIndex + 80 > indexNum >  startingIndex + 20:
		GTDic = IORobot.loadContigsFromFile(folderName, "GTDic.fasta")
		str5 = GTDic["Segkk0"]
		print startingIndex , indexNum - startingIndex
		print startingIndex, str5[startingIndex:startingIndex+100]
		print indexNum, str5[indexNum-10:indexNum+10]
		print len(segList)
		
		segList.append(str5[startingIndex:startingIndex+100])
	'''

	IORobot.writeSegOut(segList[0:lenMax], folderName, filename+ ".fasta")
	endThres = 2
	#cline = ClustalwCommandline("clustalw2", infile=folderName + filename + ".fasta", pwdnamatrix="matrix.txt",TRANSWEIGHT=0,  GAPOPEN=0, GAPEXT=0)
	cline = ClustalwCommandline("clustalw2", infile=folderName + filename + ".fasta", PWDNAMATRIX="matrix.txt",pwgapopen=0,  pwgapext=0,TRANSWEIGHT=0,GAPOPEN=0, GAPEXT=0)
	stdout, stderr = cline()

	align = AlignIO.read(folderName + filename +".aln", "clustal")
	summary_align = AlignInfo.SummaryInfo(align)
	consensus = summary_align.gap_consensus(threshold=0)

	for eachalign in align:
		if eachalign.id == 'Segkk0':
			myseq = eachalign.seq

	ctTemplate = startingIndex 
	modiList  = []
	for i in range(len(consensus)):
		#if consensus[i] == 'X':
		#	print ctTemplate, consensus[i], myseq[i], i 
		#if startingIndex + 80 > indexNum >  startingIndex + 20 and i  == 83:
		#	print  consensus[i] != myseq[i] and consensus[i] != 'X'
		#	print consensus[i], myseq[i], ctTemplate, i
		#	assert(False)
		#	print ctTemplate, consensus[i], myseq[i], i 

		if consensus[i] != myseq[i] and consensus[i] != 'X':
			if consensus[i] == '-' :
				modiList.append([ctTemplate - startingIndex, ctTemplate , 'd']) 
			elif myseq[i] == '-':	
				if len(modiList) > 1 and modiList[-1][0] == ctTemplate:
					prevIndex = int(modiList[-1][1].split('_')[-2])
					suffix =  str(prevIndex + 1)
				else: 
					suffix = "0"

				modiList.append([ctTemplate - startingIndex, ctTemplate, 'i_' + suffix+ "_"+ str(consensus[i])  ])
			else:
				modiList.append([ctTemplate - startingIndex, ctTemplate, 's_' + str(consensus[i])])
				#print consensus, i , consensus[i], myseq[i]

		if myseq[i] != '-':
			ctTemplate += 1

	modiList.sort()
	newModiList = []

		
	for k in range(len(modiList)):
		eachitem = modiList[k]
		if 20< eachitem[0]< 80:
			newModiList.append([eachitem[1], eachitem[2]])

		elif k < len(modiList) -1 and modiList[k][0] <= 20  and modiList[k+1][0] > 20 and abs(modiList[k+1][0] - modiList[k][0])  < endThres  :
			newModiList.append([eachitem[1], eachitem[2]])
		elif k >  0 and modiList[k][0] >= 80 and  modiList[k-1][0] < 80 and abs(modiList[k-1][0] - modiList[k][0])  < endThres  :
			newModiList.append([eachitem[1], eachitem[2]])

	'''
	indexNum = 10259 
	if   startingIndex + 80 > indexNum >  startingIndex + 20:
		command = "cat " + folderName + "seq.aln >> ./happy "
		os.system(command)
		print startingIndex , indexNum - startingIndex
		print consensus
		print newModiList
		#assert(False)

	'''
	return newModiList

def concatAndFormat(folderName):
	assert(True)

def computeScore(folderName, eachMatching, lambdas, interiors, readMatching, constants, isDebug, mummerLink):
	'''
	Input : lambdas, interiors, readMatching, constants
	Output : score \in real
	Algorithm : 
	1) Compute the total edits scores (from interiors, readMatching = [read2templateDic, template2readDic], and RList.fasta)
	2) Compute the abundance scores (from lambda)
	3) Combine them to give the final score 
	'''
	score = 0 
	# 1)
	editScore = 0
	q = 0.01

	if not isDebug:			

	#a) Perform an alignment and parse the results 
		alignerRobot.useMummerAlign(mummerLink, folderName, "interiorAnchor", "RList_Double.fasta", "interiors.fasta", False, "", False)
		readAnchorDic = {}
		dataList = alignerRobot.extractMumData(folderName, "interiorAnchor" +"Out")
		thres = 30
		dataList.sort(key = itemgetter(-2))
		for key, items in groupby(dataList, itemgetter(-2)):
			maxMatch = 0 
			for eachitem in items:
				if key in   readMatching[0][key] and len(readMatching[0][key]) > 0  and eachitem[-1] == readMatching[0][key][0]:
					if eachitem[4] > maxMatch :
						maxMatch = eachitem[4]
						readAnchorDic[key] = [eachitem[0], eachitem[1], eachitem[2], eachitem[3]]

		with open(folderName +  "readAnchorDic.json", 'w') as outfile:
			json.dump(readAnchorDic, outfile)
		
	else:
		readAnchorDic = readInJSON(folderName, "readAnchorDic.json")

	#b) Perform careful edit distance computation 
	interiorsDic = IORobot.loadContigsFromFile(folderName, "interiors.fasta")
	readsDic =  IORobot.loadContigsFromFile(folderName, "RList_Double.fasta")

	for i in range(len(interiors)):
		tmpScore = 0 
		for eachitem in readMatching[1]["Segkk" + str(i)]:
			readName =  eachitem[-2]
			#print readName
			if readName  in   readAnchorDic:
				readStart, readEnd, templateStart, templateEnd  = readAnchorDic[readName]
				tmpScore += Levenshtein.distance(readsDic[readName][readStart-1:readEnd] , interiorsDic["Segkk" + str(i)][templateStart-1:templateEnd])

		editScore += math.log(1.0*q/(1-2*q)) *tmpScore

	# 2)
	### Need to correct the errors

	NiList = []

	internalReads = []
	for i in range(len(lambdas)):
		internalReads += readMatching[1]["Segkk" + str(i)][0]

	internalReadsSet = set(internalReads)

	contigToReadsDic = readInJSON(folderName, "contigToReadsDic.json")
	
	for i in range(len(eachMatching)):
		leftContig, rightContig = convertName(eachMatching[i][0]), convertName(eachMatching[i][1])
		Ni = len(set(contigToReadsDic[leftContig]) - internalReadsSet) + \
			 len(set(contigToReadsDic[rightContig]) - internalReadsSet) + \
			 len(readMatching[1]["Segkk" + str(i)])
		NiList.append(Ni)


	LiList = []
	contigsDic = IORobot.loadContigsFromFile(folderName,"improved3_Double.fasta")

	for i in range(len(eachMatching)):
		leftContig, rightContig = eachMatching[i][0], eachMatching[i][1]
		left, middle, right = contigsDic[leftContig], interiorsDic["Segkk"+ str(i)], contigsDic[rightContig]
		totalLen = len(left) + len(middle) + len(right)
		overlap = IORobot.align(left, middle ,folderName,mummerLink)
		totalLen += overlap[0]
		overlap = IORobot.align(middle, right ,folderName,mummerLink)
		totalLen += overlap[0]
		LiList.append(totalLen)
	
	abunScore = 0

	for i in range(len(lambdas)):
		abunScore += math.log(lambdas[i]/LiList[i]) * NiList[i]

	# 3)
	score = editScore + abunScore
	print score, editScore, abunScore, lambdas

	return score

def convertName(myName):
	print ""
	dataInfo = myName.split("_")
	return "Segkk" + dataInfo[0][6:]

def initMatching(folderName, eachMatching, mummerLink):
	'''
	Input: eachMatching,  
	Output: prevIteration = [lambdas, templates, score] , constants = [basesMappedToEachContig]
	Algorithm: 
	1) V Init score = 0
	2) V Init templates to be concatentations 
	3) V Init lambdas to be ratio without considering the interiors
	4) V return data 
	'''	

	
	# 1) 
	score  =  0 
	
	# 2, 3) 
	templates = [] 
	basesMappedToEachContig = []
	lambdas  = []

	contigDic = IORobot.loadContigsFromFile(folderName, "improved3_Double.fasta")
	templateDic = IORobot.loadContigsFromFile(folderName, "intermediate.fasta")
	myCountDic = readInJSON(folderName, "myCountDic.json")
	lenDic = IORobot.obtainLength(folderName, "improved3_Double.fasta")

	templateName = ""
	for eachitem in templateDic:
		templateName= eachitem

	thres = 20000

	#noiselessTemplatesDic = IORobot.loadContigsFromFile(folderName, "intermediateNoiseless.fasta")
	noiselesstemplates = []

	for eachpair in eachMatching:		
		lenL, lenR = len(contigDic[eachpair[0]]) -1 , len(contigDic[eachpair[1]]) -1
		leftSeg, middleSeg, rightSeg = contigDic[eachpair[0]][-min(thres, lenL):],  templateDic[templateName],  contigDic[eachpair[1]][0:min(thres, lenR)]
		overlap1 = IORobot.align(leftSeg, middleSeg , folderName, mummerLink) 
		overlap2 = IORobot.align( middleSeg, rightSeg , folderName, mummerLink) 
		#print "overlap1 , overlap2", overlap1 , overlap2, eachpair
		#print len(rightSeg), min(thres, lenR)

		templates.append(leftSeg[0: max(len(leftSeg)-overlap1[0], 0)] + middleSeg + rightSeg[overlap2[-1]:])
		noiselesstemplates.append(leftSeg[0:max(len(leftSeg)-overlap1[0], 0)] + middleSeg + rightSeg[overlap2[-1]:])
		
		inName, outName = convertName(eachpair[0]), convertName(eachpair[1])
		numBases = np.sum(myCountDic[inName])*lenDic[eachpair[0]] + np.sum(myCountDic[outName])*lenDic[eachpair[1]]
		basesMappedToEachContig.append(numBases)


	sumOfAllBases = np.sum(basesMappedToEachContig)

	for numbase in basesMappedToEachContig:
		lambdas.append(numbase*1.0/sumOfAllBases)

	# Format return 
	IORobot.writeSegOut(noiselesstemplates, folderName, "GTDic.fasta")
	prevIteration, constants = [lambdas, templates, score] , [basesMappedToEachContig] 


	return prevIteration, constants

def formatOutputAndFinalDecision(folderName, matchingResultList):
	'''
	Input : matchingResultList = [ [prevIteration= [lambdas, interiors, score], eachMatching] , ... ]
	Output : score, matching, contentForBetterInteriorToFlank
	Algorithm : 
	1) Select the one with largest score
	2) Format and return 
	Heuristic here for scoring though... should improve later

	'''
	score, matching, contentForBetterInteriorToFlank = 0, 0 ,0 
	scoreList = []


	maxScore = -10**9
	for eachitem in matchingResultList:
		print "eachMatching, lambdas, scores", eachitem[1], eachitem[0][0], eachitem[0][-1]
		eachMatching, lambdas, scores = eachitem[1], eachitem[0][0], eachitem[0][-1]
		scoreList.append(scores)

		if scores > maxScore: 
			score, matching, contentForBetterInteriorToFlank = eachitem[0][-1], eachitem[1],  eachitem[0][1]

	scoreList.sort(reverse=True)

	ratioScore = scoreList[0]/scoreList[1]

	return ratioScore, matching, contentForBetterInteriorToFlank

def formMathchingList(folderName, contigLeftFileName, contigRightFileName):
	matchingList = []

	LNameList = readInJSON(folderName, contigLeftFileName)
	RNameList = readInJSON(folderName, contigRightFileName)
	
	for eachitem in  permutations(RNameList):
		RList = list(eachitem)
		tmpList = []
		for i in range(len(RList)):
			tmpList.append([LNameList[i], RList[i]])

		matchingList.append(tmpList)

	return matchingList
	
def generateAssociatedReadDic(folderName):
	dataList = []
	numberOfFiles = houseKeeper.globalParallelFileNum
	for i in range(1, 1+numberOfFiles): 
		if i < 10:
			indexOfMum = "0" + str(i)
		else:
			indexOfMum = str(i)
		
		dataList = dataList+ alignerRobot.extractMumData(folderName, "outAbun"+ str(indexOfMum)+"Out")

	dataList.sort(key=itemgetter(-1))

	contigToReadsDic = {}

	lenContigDic = IORobot.obtainLength(folderName , "improved3.fasta")
	for eachitem in lenContigDic:
		contigToReadsDic[eachitem] = []

	for key, items in groupby(dataList, itemgetter(-1)):
		maxLen = 0
		tmpTarget = ""
		for eachitem in items:
			if eachitem[-4] > maxLen : 
				maxLen = eachitem[-4]
				tmpTarget = eachitem[-2]

		contigToReadsDic[tmpTarget].append(key)

	with open(folderName +  "contigToReadsDic.json", 'w') as outfile:
		json.dump(contigToReadsDic, outfile)
	
def debugging(folderName):
	referenceDic = IORobot.loadContigsFromFile(folderName, "reference.fasta")
	interiorsDic = IORobot.loadContigsFromFile(folderName, "interiors.fasta")
	GTDic = IORobot.loadContigsFromFile(folderName, "GTDic.fasta")

	str1 = referenceDic["Segkk0"][2500000:2500000+12000]
	str2 = referenceDic["Segkk1"][2500000:2500000+12000]
	print Levenshtein.distance(str1, str2)
	print Levenshtein.editops(str1, str2)

	str3 = interiorsDic["Segkk0"][7000:7000+12000]
	str4 = interiorsDic["Segkk1"][7000:7000+12000]
	print Levenshtein.distance(str1, str4)
	print Levenshtein.distance(str2, str3)
	print Levenshtein.editops(str1, str4)
	
	print ""
	offset = 4000
	print Levenshtein.editops(str2, str3)
	print str1[offset-10:offset+10]
	print str2[offset-10:offset+10]
	print str3[offset-10:offset+10]
	print str4[offset-10:offset+10]

	str5 = GTDic["Segkk0"][7000:7000+12000]
	str6 = GTDic["Segkk1"][7000:7000+12000]
	
	print str5[offset-10:offset+10]
	print str6[offset-10:offset+10]
			

	print Levenshtein.editops(str2, str4)

def loadRListDic(folderName):
	numberOfFiles = houseKeeper.globalParallelFileNum
	thres = 10000

	dataList = []
	for i in range(1, 1+numberOfFiles): 
		if i < 10:
			indexOfMum = "0" + str(i)
		else:
			indexOfMum = str(i)

		dataList = dataList+ alignerRobot.extractMumData(folderName, "outAbun"+ str(indexOfMum)+"Out")

	dataList.sort(key = itemgetter(-2))

	lenDic = IORobot.obtainLength(folderName, "improved3.fasta")

	RListDic = {}
	for key, items in groupby(dataList, itemgetter(-2)): 
		RListDic[key] = []
		for eachitem in items:
			if eachitem[0] <  thres or eachitem[1] >lenDic[eachitem[-2]] - thres: 
				readName = eachitem[-1]
				#RListDic[key].append("Contig" + readName[5:] + "_p")
				#RListDic[key].append("Contig" + readName[5:] + "_d")
				RListDic[key].append(readName)
		RListDic[key] = abunHouseKeeper.getDistinct(RListDic[key])

	return RListDic

def XResolvePreparation(Gnew, GContigRead, Grev, folderName, myCountDic, lenDic, N1, mummerLink):
	print "XResolvePreparation"
	'''
	Implementation steps (TODO : @ kakitfive): 
		1) From the log files, generate :
			contigLeft.json [automatic from XResolve] [Issue : what if there is mulitple to aggregate?] V
			contigRight.json [automatic from XResolve] [Issue : what if there is mulitple to aggregate?] V
			RList.fasta	[in preparation() and generateAssociatedReadDic()] V 
			intermediate.fasta [The most tricky part]
				a) Find multiple paths u->X->v (s.t. it is a robust one ?! Need to separate into cases for XResolve and BResolve)
				b) Find common region as the intermediate.fasta (is it good enough ?! )
					i) Will it work in the synthetic case? 
					ii) Will it work in the real case with chimeras and other noisy parts ? 
					iii) Is it the simplest method that you can think of ? 
					iv) Draw a conclusion and note the strengths and weaknesses

		2) Run the EMFlow

		3) Write down the results

		4) Rewire the interface to make it adaptable

		5) Passs the interiors back as well

	'''
	# output format : combinedList [[], [], [[0, 8]], [[5, 7]], [], [], [], [], [], []]
	RListDic = loadRListDic(folderName)

	xResolvedList = [ [] for i in range(N1)]

	contigDic = IORobot.loadContigsFromFile(folderName, "improved3_Double.fasta")

	for nodeI in range(N1):
		v = Gnew.graphNodesList[nodeI]

		inList = [] 
		for eachitem in v.listOfPrevNodes:
			inList.append(abunHouseKeeper.parseIDToName(eachitem[0],'C', 0))

		outList = []
		for eachitem in v.listOfNextNodes:
			outList.append(abunHouseKeeper.parseIDToName(eachitem[0],'C', 0))

		print "inList, outList", inList, outList

		if len(inList) > 1 and len(outList) > 1:
			with open(folderName +  "contigLeft.json", 'w') as outfile:
				json.dump(inList, outfile)
			
			with open(folderName +  "contigRight.json", 'w') as outfile:
				json.dump(outList, outfile)

			#xnodeName = "Segkk" + str(nodeI/2)
			xnodeName = "Segkk" + str(nodeI/2)

			IORobot.writeSegOut([contigDic[abunHouseKeeper.parseIDToName(nodeI,'C',0)]], folderName, "intermediate.fasta")

			RList = RListDic[xnodeName]
			IORobot.putListToFileO(folderName, "raw_reads.fasta", "RList", RList)
			IORobot.writeToFile_Double1(folderName, "RList.fasta", "RList_Double.fasta", "contig")

			score, matching, contentForBetterInteriorToFlank = EMFlow(folderName, mummerLink)
			ansList = []

			if 1/score > 1.001 : 
				for eachsub in matching:
					ansList.append([abunHouseKeeper.parseEdgeNameToID(eachsub[0], 'C'), abunHouseKeeper.parseEdgeNameToID(eachsub[1], 'C')])

			xResolvedList[nodeI] = ansList

	return xResolvedList

def findPathList(folderName, G, N1, inList, outList):
	# return the list of paths
	pathList = [] 
	inIndexList, outIndexList = [], []

	for i in inList:
		inIndexList.append(abunHouseKeeper.parseEdgeNameToID(i,'C'))
	
	for i in outList:
		outIndexList.append(abunHouseKeeper.parseEdgeNameToID(i,'C'))
	

	for leftCtgIndex in inIndexList:
		for rightCtgIndex in outIndexList:
			#tmpPathsList = abunGraphLib.findPathBtwEndsFast(folderName, leftCtgIndex, rightCtgIndex, G, N1)
			tmpPathsList = abunGraphLib.BFS(leftCtgIndex, rightCtgIndex, G, N1)
			#print tmpPathsList
			#assert(False)
			pathList.append( tmpPathsList)

	return pathList 

def findAPair(pathList):
	# find a mutually exclusive path that start/end differently
	#print pathList
	#assert(False)
	paths = [] 
	ck = False 
	for i in range(len(pathList) - 1):
		for j in range(i, len(pathList)):
			#print i, j, pathList[i], pathList[j]
			if pathList[i][0] != pathList[j][0] and pathList[i][-1] != pathList[j][-1]:
				paths = [[pathList[i][0],pathList[i][-1] ] , [pathList[j][0], pathList[j][-1]]]
				ck = True
				break
		if ck: 
			break
	#print paths
	#assert(False)
	return paths 

def findPathSegments(folderName, paths, N1, mummerLink):
	# find the segments by overlap and joining subroutine in IORobot
	contigFilename = "improved3"
	readsetFilename = "phasingSeedName"
	contigReadGraph = "phaseStringGraph1"

	contigsDic = IORobot.loadContigsFromFile(folderName, "improved3_Double.fasta")

	pathSegList = [] 
	for eachpath in paths:
		tmpSeq = ""
		leftCtgIndex ,rightCtgIndex, leftEnd, rightStart, middleContent = \
			abunGraphLib.singleGapLookUp(eachpath ,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,"phasingSeedName")
		print "Now",  leftCtgIndex, rightCtgIndex, eachpath,leftEnd, rightStart, len(middleContent)

		#assert(False)

		tmpSeq = contigsDic[abunHouseKeeper.parseIDToName(leftCtgIndex,'C', 0)][max(0, leftEnd-20000):leftEnd] + \
					middleContent + contigsDic[abunHouseKeeper.parseIDToName(rightCtgIndex,'C', 0)][rightStart:min(rightStart + 20000, len(contigsDic[abunHouseKeeper.parseIDToName(rightCtgIndex,'C', 0)]) -1)]

		pathSegList.append(tmpSeq)


	return pathSegList[0], pathSegList[1]

def mainFlow():
	'''
	Input : RList.fasta, contigLeft.json, contigRight.json, intermediate.fasta
	Output : score, matching, betterInteriorToFlank.fasta 
	'''
	print "Perform EM Algorithm"
	#folderName = "dataFolder/"
	
	folderName = "brtest/"

	if False:
		XResolvePreparation(folderName)

	if False:
		BResolvePreparation(folderName)

	if False:
		generateAssociatedReadDic(folderName) 	
		preparation(folderName)

	if False:
		EMFlow(folderName)

	if False:
		debugging(folderName)

	if True:
		contigReadGraph = "phaseStringGraph1"
		G = abunGraphLib.seqGraphWt(0)
		G.loadFromFile(folderName, contigReadGraph)
		N1 = 8
		Grev = abunGraphLib.formReverseGraphFast(G)

		#inList , outList = [0, 8] , [5, 13] 
		inList , outList =  [6, 14], [3, 11] 
		BResolvePreparation(folderName, inList, outList, G, Grev, N1)

def EMFlow(folderName, mummerLink):
	matchingList = [] 
	numberOfIterations = 1
	isDebug = False
	matchingResultList = []
	
	matchingList = formMathchingList(folderName, "contigLeft.json", "contigRight.json")
	matchingList.sort(reverse = True)
	# print matchingList

	for eachMatching in matchingList:
		prevIteration, constants = initMatching(folderName, eachMatching, mummerLink)

		for iterationI in range(numberOfIterations):	
			readMatching = computeReadAssociation(folderName, prevIteration, constants, isDebug, mummerLink)
			lambdas = computeLambdaList(folderName, prevIteration, readMatching, constants)
			interiors = computeInteriorList(folderName , prevIteration, readMatching, constants, isDebug, mummerLink )
			score = computeScore(folderName, eachMatching, lambdas, interiors, readMatching, constants, isDebug, mummerLink)
			prevIteration = [lambdas, interiors, score] 

		#assert(False)
		matchingResultList.append([prevIteration, eachMatching])

	score, matching, contentForBetterInteriorToFlank = formatOutputAndFinalDecision(folderName, matchingResultList)

	return score, matching, contentForBetterInteriorToFlank

def	BResolvePreparation(folderName, inList, outList, G, Grev, N1, mummerLink):
	print "BResolvePreparation"
	# format :  resolvedList, brResolvedList, inList, outList [] [[3, 1], [3, 7]] [6] [3, 15]
	# resolvedList in standard format ... just that inList, outList has unnecessary *2 for head/tail difference
	# Input : brtest/ [0, 8] [5, 13]
	# print folderName,  inList, outList

	resolvedList = []
	
	if len(inList) > 1 and len(outList) > 1:
		# prepare left/righ contigs 
		contigLeft = []
		for eachitem in inList:
			contigLeft.append(abunHouseKeeper.parseIDToName(eachitem/2,'C', 0))

		contigRight = []
		for eachitem in outList:
			contigRight.append(abunHouseKeeper.parseIDToName(eachitem/2,'C', 0))

		print "contigLeft, contigRight", contigLeft, contigRight

		with open(folderName +  "contigLeft.json", 'w') as outfile:
			json.dump(contigLeft, outfile)

		with open(folderName +  "contigRight.json", 'w') as outfile:
			json.dump(contigRight, outfile)

		# prepare RList
		RListDic = loadRListDic(folderName)

		RList = [] 
		for eachkey in contigLeft + contigRight:
			nodeIndex = abunHouseKeeper.parseEdgeNameToID(eachkey, 'C')
			nodeName = "Segkk" + str(nodeIndex/2)
			RList = RList + RListDic[nodeName]

		RList  = abunHouseKeeper.getDistinct(RList)

		IORobot.putListToFileO(folderName, "raw_reads.fasta", "RList", RList)
		IORobot.writeToFile_Double1(folderName, "RList.fasta", "RList_Double.fasta", "contig")

		# prepare intermediate
		### Look for a path and then join here. 
		
		pathList = findPathList(folderName, G, N1, contigLeft, contigRight)
		paths = findAPair(pathList)
		path1, path2 = findPathSegments(folderName, paths, N1, mummerLink)

		IORobot.writeSegOut([path1], folderName, "path1.fasta")
		IORobot.writeSegOut([path2], folderName, "path2.fasta")

		alignerRobot.useMummerAlign(mummerLink, folderName, "comparison", "path1.fasta", "path2.fasta", False, "", False)
		dataList = alignerRobot.extractMumData(folderName, "comparison" +"Out")

		dataList.sort(key=itemgetter(-2))
		begin, end = 1 , 100

		for key, items in groupby(dataList, itemgetter(-2)):
			maxLen = -1 
			for eachitem in items:
				if eachitem[4] > maxLen:
					begin, end = eachitem[0], eachitem[1] 
					maxLen = eachitem[4]

		path1Dic = IORobot.loadContigsFromFile(folderName, "path1.fasta")

		IORobot.writeSegOut([path1Dic["Segkk0"][begin-1:end]], folderName, "intermediate.fasta")

		ratioScore, matching, contentForBetterInteriorToFlank = EMFlow(folderName, mummerLink)

		#assert(False)
		if 1/ratioScore > 1.001:
			print "kkbug score", ratioScore 
			for eachsub in matching:
				resolvedList.append([abunHouseKeeper.parseEdgeNameToID(eachsub[0], 'C'), abunHouseKeeper.parseEdgeNameToID(eachsub[1], 'C')])

	return resolvedList


#t0 = time.time()
#mainFlow()
#print "Time spent (s) : ", time.time() - t0

