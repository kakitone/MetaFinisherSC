import matplotlib.pyplot as plt
import json

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import nonRedundantResolver

from itertools import groupby
from operator import itemgetter

import abunGraphLib
from finisherSCCoreLib import graphLib
import abunHouseKeeper
import abunSplitter

import argparse
import os

import time
import datetime
import re


def test1():
	lenDic = {}
	coverageDic = {}
	
	lenDic = IORobot.obtainLength("/Users/kakitlam/", "abun.fasta")
	
	f = open("/Users/kakitlam/Documents/abundata", 'r')
	tmp = f.readline()
	
	while len(tmp) > 0:
		if len(tmp) > 10:
			myitem = tmp[0:-1].split()
			coverageDic[myitem[0]] = float(myitem[1])
		tmp = f.readline()
	
	f.close()
	
	myList = []
	baseCt = {}
	
	for eachitem in lenDic:
		myList.append(lenDic[eachitem]*coverageDic[eachitem])
		baseCt[eachitem] = lenDic[eachitem]*coverageDic[eachitem]
	
	
	for eachitem in lenDic :
		print eachitem,  baseCt[eachitem]
	
	
	
	for eachitem in lenDic :
		print eachitem, lenDic[eachitem]
	
	
	for eachitem in lenDic :
	        print eachitem, coverageDic[eachitem]


def test2():
	import abunSplitter
	abunSplitter.mainFlow("Apr10Test/", "/usr/bin/")

def viewLenDic():
	
	folderName = "Apr10Test/"
	json_data = open(folderName + "myCountDic.json", 'r')
	myCountDic = json.load(json_data)
	
	contigLenDic = IORobot.obtainLength(folderName,  "LC_n.fasta")
	
	toPlotListX = []
	toPlotListY = []
	
	for eachitem in contigLenDic:
		toPlotListX.append(myCountDic[eachitem])
		toPlotListY.append(contigLenDic[eachitem])
	
	print toPlotListX, toPlotListY
	
	
	with open(folderName + "toPlotListX.json", 'w') as f:
		json.dump(toPlotListX, f)
		
	
	with open(folderName + "toPlotListY.json", 'w') as f:
		json.dump(toPlotListY, f)
	#plt.scatter(toPlotListX, toPlotListY)
	
	#plt.savefig(folderName + 'foo.png')	
		
		
def coloringNodes():
	folderName = "Apr10Test/"
	if False:
		alignerRobot.useMummerAlign("/usr/bin/", folderName, "debug", "reference.fasta", "LC_n.fasta")
	
	dataList = alignerRobot.extractMumData(folderName, "debugOut")
	
	dataList.sort(key = itemgetter(-1))
	
	mappedDic = {}
	
	for key, items in groupby(dataList, itemgetter(-1)):
		print "key", key
		matchLen = -1
		
		for eachitem in items: 
			if eachitem[-4] > matchLen:
				mappedDic[key]  = eachitem[-2]
				matchLen = eachitem[-4]

	
	for eachitem in mappedDic:
		if mappedDic[eachitem] == 'c3':
			print str(int(eachitem[5:])*2)+"{color:blue}"
			print str(int(eachitem[5:])*2+1)+"{color:blue}"
		
		if mappedDic[eachitem] == 'c1':
			print str(int(eachitem[5:])*2)+"{color:green}"
			print str(int(eachitem[5:])*2+1)+"{color:green}"

def coloringNodes2():
	folderName = "Apr10Test/"
	
	json_data = open(folderName + "myCountDic.json", 'r')
	myCountDic = json.load(json_data)
	
	for eachitem in myCountDic:
		if myCountDic[eachitem] > 200:
			print str(int(eachitem[5:])*2)+"{color:blue}"
			print str(int(eachitem[5:])*2+1)+"{color:blue}"
		
		if myCountDic[eachitem] <= 200:
			print str(int(eachitem[5:])*2)+"{color:green}"
			print str(int(eachitem[5:])*2+1)+"{color:green}"
	

def mapStrangePairs():
	folderName = "Apr10Test/"
	
	json_data = open(folderName + "furtherGapList.json", 'r')
	furtherGapList = json.load(json_data)
	
	segLookUp = IORobot.readContigsFromFile(folderName, "LC_n_Double.fasta")
	
	f = open(folderName + "wrongCondense.fasta", 'w')
	ctr = 0
	for eachitem in furtherGapList:
		beforeI, afterI = eachitem[0], eachitem[1]
		
		f.write(">Segkk"+str(ctr)+"\n")
		f.write(segLookUp[beforeI]+"\n")
		ctr = ctr + 1 
		
		f.write(">Segkk"+str(ctr)+"\n")
		f.write(segLookUp[afterI]+"\n")
		ctr = ctr + 1 
	
	f.close()
	
	if False:
		alignerRobot.useMummerAlign("/usr/bin/", folderName, "wrongCondenseDebug", "reference.fasta", "wrongCondense.fasta")
	
	dataList = alignerRobot.extractMumData(folderName, "wrongCondenseDebugOut")
	
	dataList.sort(key = itemgetter(-1))
	
	mappedDic = {}
	
	for key, items in groupby(dataList, itemgetter(-1)):
		print "key", key
		matchLen = -1
		
		for eachitem in items: 
			if eachitem[-4] > matchLen:
				mappedDic[key]  = eachitem
				matchLen = eachitem[-4]
				
	
	for eachitem in mappedDic:
		print "results : ", eachitem, mappedDic[eachitem]


def continuousIntegration():
	if False:
		G = graphLib.seqGraph(10)
		for i in range(5):
			G.insertEdge(i,i+1,1997)
			G.insertEdge(i,i+2, 1997)

		resultList = abunGraphLib.BFS_revisit(1,3,G,1)

		print "resultList", resultList 

	if False : 

		folderName, mummerPath, directPathList, indirectPathList, contigFile, readFile = \
			"Apr10Test/", "/usr/bin/", [[1, 486, 217], [1, 8642, 217], [1, 13465, 217]], [[1, 486, 217]], "improved3_Double.fasta", "phasingSeedName_Double.fasta"

		abunGraphLib.formPathSeq(folderName, mummerPath, directPathList, indirectPathList, contigFile, readFile)
    
		if False:
			lenDic = IORobot.obtainLength(folderName , contigFile)
			N1 = len(lenDic)

			print "N1", N1

			G = graphLib.seqGraph(0)
			G.loadFromFile(folderName, "phaseStringGraph1")

			adj = [[] for i in range(N1)]

			for i in range(N1): 
			    adj[i] = abunGraphLib.findAllReachable(i, N1, G)

			Gnew = abunGraphLib.seqGraphDynamic(N1)

			for i in range(N1):
			    for j in adj[i]:
			        Gnew.insertEdge(i,j,1997)


			Gnew.initAdv()    
			Gnew.doubleEdgeReduction()

			contigPaths = abunGraphLib.findAllPathK(1, 217, Gnew, 3)
			contigReadPaths = abunGraphLib.findAllPathK(1, 217, G, 5)

			print "contigPaths", contigPaths
			print "contigReadPaths", contigReadPaths

			Gnew.transitiveReduction()

	if False:
		toDelete = abunGraphLib.decideCut("Apr10Test/", "/usr/bin/")
		print toDelete

	if False:
		G = graphLib.seqGraph(0)
		G.loadFromFile("Apr10TestA/", "xResolvedGraph")

		if False:
			for i in range(len(G.graphNodesList)):

				v = G.graphNodesList[i]

				if len(v.nodeIndexList) > 0:
					print i , v.listOfPrevNodes , v.listOfNextNodes

		G.reportEdge()
		lenDic = IORobot.obtainLength("Apr10TestA/", "improved3_Double.fasta")
		mylist = [401, 207, 405, 407, 344]

		json_data = open("Apr10TestA/" + "myCountDic.json", 'r')
		myCountDic = json.load(json_data)

		for x in mylist:
			print x, lenDic["Contig"+str(x/2)+"_p"], myCountDic["Segkk"+str(x/2)]


	if False:
		folderName = "Apr10TestA/"
		G = graphLib.seqGraph(0)
		G.loadFromFile(folderName , "xResolvedGraph")

		json_data = open(folderName + "mapDummyToRealDic.json", 'r')
		mapDummyToRealDic = json.load(json_data)

		lenDic = IORobot.obtainLength(folderName, "improved3_Double.fasta")
		print len(G.graphNodesList)
		print len(mapDummyToRealDic)
		
		print "fake N1 , real N1 ", len(G.graphNodesList) - len(mapDummyToRealDic), len(lenDic)


	if False:
		abunSplitter.mainFlow("Apr10TestB/", "/usr/bin/")

	if False: 
		nonRedundantResolver.removeEmbedded("Apr10TestD/", "/usr/bin/")

	if False:
		folderName, contigReadGraph = "Apr10TestA/", "phaseStringGraph1"
		G = graphLib.seqGraph(0)
		kthres, edgeThres = 3, 1
		G.loadFromFile(folderName, contigReadGraph)
		lenDic = IORobot.obtainLength(folderName , "improved3_Double.fasta")

		N1 = len(lenDic)

		adj = [[] for i in range(N1)]

		for i in range(N1): 
		    tmpList = abunGraphLib.findAllReachable(i, N1, G)
		    
		    for j in tmpList:
		        if len(abunGraphLib.findAllPathK(i,j,G,kthres)) >= edgeThres:
		            adj[i].append(j) 

		    #print i, adj[i]

	    ### Filter adaptor skipped case 

		adaptorPair = []

		for i in range(len(adj)):
		    if  i % 2 == 0:
		        if i + 1 in adj[i]:
		            adj[i].remove(i+1)
		            adaptorPair.append([i, i+1])
		    elif i % 2 ==1: 
		        if i-1 in adj[i] :
		            adj[i].remove(i-1)
		            adaptorPair.append([i, i-1])

		Gnew = abunGraphLib.seqGraphDynamic(N1)

		for i in range(N1):
		    for j in adj[i]:
		        Gnew.insertEdge(i,j,1997)

		for eachpair in adaptorPair:
		    u, v = eachpair[0], eachpair[1]
		    for x in Gnew.graphNodesList[u].listOfPrevNodes:
		        xIndex = x[0]
		        Gnew.removeEdge(xIndex, v)
		    for y in Gnew.graphNodesList[v].listOfNextNodes:
		        yIndex = y[0]
		        Gnew.removeEdge(u, yIndex)


        #Gnew.reportEdge()
		count2 = 0
		for i in range(len(Gnew.graphNodesList)):
			if  len(Gnew.graphNodesList[i].listOfPrevNodes) == 2 and  len(Gnew.graphNodesList[i].listOfNextNodes) == 2:
				count2 = count2 + 1
				print str(i)+"{color:red}"

		print "count2, ", count2

		### End filter adaptor skipped case 
	if True:
		nonRedundantResolver.removeRedundantWithFile("May11TestB/" , "/usr/bin/", "abun", "abunDebug", "abunNoEmbed")





def multipleTesting(folderName, mummerLink, option):
	print "Experiments : "
		
	header = "/data/kakitone/May07-2015/workingMetaFinisherSC/"

	os.system("mkdir "+ header + folderName )
	
	os.system("cp "+header+"Apr10Test/*  "  + header + folderName )

	if option == "1":
		abunHouseKeeper.abunGlobalSplitParameterRobot.runGraphSurgery = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runXResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.RThres = 5
		abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresB  = 2
		#abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresX  = 1
		abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerB = 1.5
		abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperB = 1.5
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerX = 0.5
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperX = 1.5

	elif option == "2" : 
		abunHouseKeeper.abunGlobalSplitParameterRobot.runGraphSurgery = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runXResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.RThres = 5
		abunHouseKeeper.abunGlobalSplitParameterRobot.BRThres  = 2
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerB = 2
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperB = 2
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerX = 0.5
		#abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperX = 0.5
	elif option == "3": 

		abunHouseKeeper.abunGlobalSplitParameterRobot.runGraphSurgery = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.runXResolve = True
		abunHouseKeeper.abunGlobalSplitParameterRobot.RThres = 5
		abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresB = 3
		abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresX = 1



	abunSplitter.mainFlow(folderName, mummerLink)

	os.system("python /data/kakitone/download2/quast-2.3/quast.py " + header + folderName +"abun.fasta -o " + header + folderName + "   -R  "+ header + folderName + "reference.fasta")

def testingsymmetry():
	print "Hello world"

def loggHeaders():
	folderName, filename= "/data/kakitone/May07-2015/workingMetaFinisherSC/", "results.csv"

	myDic =  abunHouseKeeper.abunGlobalSplitParameterRobot.__dict__
	f = open(folderName + filename , 'a')
	tmpString  = ""
	st = "time"
	
	tmpString = tmpString + str(st) + ","

	for eachitem in myDic:
		tmpString = tmpString + eachitem + ","

	tmpString = tmpString  + "misassembly" + "," + "ncontig"+ "\n"

	f.write(tmpString)
	f.close()


def loggResults():
	folderName, filename= "/data/kakitone/May07-2015/workingMetaFinisherSC/", "results.csv"
	print "Logging results "

	myDic =  abunHouseKeeper.abunGlobalSplitParameterRobot.__dict__
	f = open(folderName + filename , 'a')
	tmpString  = ""
	st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d~%H:%M:%S')
	
	tmpString = tmpString + str(st) + ","

	for eachitem in myDic:
		tmpString = tmpString + str(myDic[eachitem]) + ","
	
	# Use regex to get misassemby and Ncontig here
	quastReport = folderName + "quast_results/latest/report.txt"
	
	regex = re.compile("misassemblies*")

	misassembly = str(1997)
	ncontig = str(1997) 
	
	with open(quastReport) as f2:
	    for line in f2:
	        result = re.match("# misassemblies (.*?) .*", line)
	        if result != None: 
		        ans = result.group(0).split()
		        misassembly = ans[2]
	        
	        result = re.match("# contigs \(>= 0 bp\) (.*?) .*", line)
	        if result != None: 
		        ans = result.group(0).split()
		        ncontig=ans[-1]

	
	tmpString = tmpString  + misassembly + "," + ncontig+ "\n"

	f.write(tmpString)

	f.close()

	
#mapStrangePairs()

#coloringNodes2()
#viewLenDic()
	
#test2()

#mapStrangePairs()

continuousIntegration() 

#loggResults()

if False:
	parser = argparse.ArgumentParser(description='aSplitter')

	parser.add_argument('folderName')
	parser.add_argument('mummerLink')
	parser.add_argument('option')

	args = vars(parser.parse_args())

	multipleTesting(args['folderName'], args['mummerLink'], args['option'])


