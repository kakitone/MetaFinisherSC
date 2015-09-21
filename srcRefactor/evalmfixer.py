
import json
import bisect 
from srcRefactor.repeatPhaserLib.finisherSCCoreLib \
	import IORobot
from srcRefactor.repeatPhaserLib.finisherSCCoreLib \
	import alignerRobot
from operator import itemgetter
from itertools import groupby
import os ,sys
import argparse
import srcRefactor.misassemblyFixerLib.intervalunion as intervalunion


'''
Goal: Calculate the precision/recall for MFixer

Input : contigs.fasta, reference.fasta 
Output: [precision, recall] at each step of MFixer

Algorithm: 
	1) GTFinder { I: LC.fasta, LC_filtered.fasta, reference.fasta; O: GTMap = [ [ contigsName,  [ [start1,end1], [start2, end2], ... ]] , ... , ] }
	2) BKPtMarker {I: intermediate files of MFixer ; O: BkPtInfo = [ [contigName, [bp1, bp2, ...,] ], ..., ]}
	3) CalculatePreRecall {I:  GTMap, BkPtInfo ; O: precision/recall}
'''

def readInJSON(folderName, filename):
	json_data = open(folderName + filename, 'r')
	dataItem = json.load(json_data)
	return dataItem

def GTFinder(folderName,inputfile):
	# GTFinder { I: LC.fasta, LC_filtered.fasta, reference.fasta; O: GTMap = [ [ contigsName,  [ [start1,end1], [start2, end2], ... ]] , ... , ] }
	# "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  |        scf7180000000702    ref_NC_001133_"

	### Finding the alignment 
	if False:
		alignerRobot.useMummerAlign("/usr/bin/", folderName, "groundTruthMatchFixer"+inputfile, inputfile, "reference.fasta", False, "", False)
	
	dataList = alignerRobot.extractMumData(folderName, "groundTruthMatchFixer" +inputfile+ "Out")

	lenDic = IORobot.obtainLength(folderName, inputfile)
	#print len(dataList)
	### Parsing the alignment 
	GTMap = []
	dataList.sort(key=itemgetter(-2))

	for key, items in groupby(dataList, itemgetter(-2)):
		contigName = key 
		tmpList = list(items)
		tmpList.sort(key = itemgetter(0))
		#print len(tmpList)
		#rangeList= rangeParser(tmpList, lenDic[key])
		thres = 100
		B = intervalunion.intervalCover(tmpList, thres)
		rangeList = intervalunion.reportMisAssemblyIntervals(B, lenDic[key], thres)

		GTMap.append([contigName, rangeList])


	with open(folderName + inputfile+ "GTMap.json", 'w') as outfile:
		json.dump(GTMap, outfile)

def rangeParser(tmpList, contigLen):
	thres  = 30
	x = thres
	rangeList = [] 
	lastStart = 0  
	
	while x <= contigLen - thres: 
		#print "x", x
		y, z, lastStart = findLRBracket(x, tmpList, lastStart, thres)
		if x == thres:
			pass
			#if y>x:
			#	rangeList.append([x,y])
		else:
			if y < x - thres:
				rangeList.append([y,x])

		if x == z:
			#print "x, y, z, contigLen : ", x, y, z, contigLen, tmpList
			#assert(x!=z)
			rangeList.append([max(1,x-thres), min(x+thres, contigLen-1)])
			break

		x = z


	return rangeList

def findLRBracket(x,dataList, lastStart, thres):
	y, z, lastStart  = 0 , 0 , 0 

	index = bisect.bisect(dataList,[x + thres] + [0 for i in range(len(dataList[0]) -1)])
	
	#print dataList[0], x+ thres, lastStart, index , len(dataList)

	endindex = min(index + 1, len(dataList))
	for i in range(lastStart, endindex):
		if dataList[i][1] > z : 
			z = dataList[i][1]
			y = dataList[i][0]

	lastStart = endindex

	return y,z,lastStart

def BKPtMarker(folderName, filename):
	# BKPtMarker {I: intermediate files of MFixer ; O: BkPtInfo = [ [contigName, [bp1, bp2, ...,] ], ..., ]}
	
	bkptInfoList = readInJSON(folderName, filename)
	#bkptInfoList = readInJSON(folderName, "repeatDic.json")
	#bkptInfoList = readInJSON(folderName, "modifiedOutliners.json")
	
	BkPtInfo = []

	for key in bkptInfoList:
		BkPtInfo.append([key, bkptInfoList[key]])
	
	with open(folderName + "BkPtInfo.json", 'w') as outfile:
		json.dump(BkPtInfo, outfile)

def CalculatePreRecall(folderName, filename, bkPtList, carryover=0):
	indComponentTestListMFixer = []
	
	for i in range(len(bkPtList)):
		GTMap = readInJSON(folderName,filename+ "GTMap.json")
		BKPtMarker(folderName, bkPtList[i])
		BkPtInfo = readInJSON(folderName, "BkPtInfo.json")

		T = carryover
		P = 0
		TP = 0
		FP = 0

		GTDic = {}

		for eachitem in BkPtInfo:
			P = P  + len(eachitem[1]) - 2

		for eachitem in GTMap:
			T = T + len(eachitem[1])
			GTDic[eachitem[0]] = eachitem[1]
			#print eachitem[0], eachitem[1]

		for eachitem in BkPtInfo:
			#print eachitem
			contigName  = eachitem[0]
			for eachbkpt in eachitem[1][1:-1]:
				
				ck = False
				if contigName in GTDic:
					for ckrange in GTDic[contigName]:
						if ckrange[0] <= eachbkpt <= ckrange[1]:
							ck = True
							GTDic[contigName].remove(ckrange)
							break


				if ck:	
					TP += 1 
				if not ck:
					FP +=1

		if P > 0 :
			precision = TP*1.0/P
		else:
			precision = 1

		recall = TP*1.0/T

		if recall != 0 :
			f1score = recall*precision*2/(precision + recall)
		else:
			f1score = 0

		print "i, precision, recall, TP_num, FP_num  : %d \t %f \t %f \t %d \t %d \t %f "  % (i, precision, recall, TP ,FP, f1score)  
		indComponentTestListMFixer.append([precision, recall, TP, FP, bkPtList[i] ])


	return indComponentTestListMFixer
	
def mainFlow(folderName):

	print "MFixer analysis."
	#folderName = os.path.abspath(os.path.dirname(sys.argv[0])) + "/" +folderName

	bkPtList = ["adaptorSkippedLogDic.json"]
	GTFinder(folderName,"LC.fasta")
	indComponentTestListMFixer = CalculatePreRecall(folderName,"LC.fasta",bkPtList, 0)


	GTFinder(folderName,"LC_filtered.fasta")
	bkPtList = ["blkDic.json",  "repeatDic.json", "modifiedOutliners.json" , "blkDicNew.json"] 
	indComponentTestListMFixer2 = CalculatePreRecall(folderName, "LC_filtered.fasta", bkPtList, indComponentTestListMFixer[0][-2])

	combined = indComponentTestListMFixer + indComponentTestListMFixer2

	with open(folderName + 'indComponentTestListMFixer.json', 'w') as f:
		json.dump(combined, f)
	
def unitTesting():
	print "unittesting"
	folderName = "dataFolder/"
	if True:
		GTFinder(folderName, "LC.fasta")
	if False:
		BKPtMarker(folderName)
	if False:
		CalculatePreRecall(folderName)

parser = argparse.ArgumentParser(description='evalmfixer')
parser.add_argument('folderName')
args = vars(parser.parse_args())

folderName = args['folderName']
mainFlow(folderName)













