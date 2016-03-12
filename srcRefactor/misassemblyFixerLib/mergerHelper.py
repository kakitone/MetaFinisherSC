from ..repeatPhaserLib.finisherSCCoreLib import IORobot
from ..repeatPhaserLib.finisherSCCoreLib import alignerRobot
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper 
import bisect
import adaptorFix
import groupctest


from operator import itemgetter
from itertools import groupby
import json
import numpy as np
import argparse
import  intervalunion
import breakPointFinding


def repeatFinder(folderName, inputName):
	'''
	Input : dataList
	Output : repeatIntervalDic
	'''
	dataList = alignerRobot.extractMumData(folderName, "self"+inputName+"Out")
	dataList = alignerRobot.transformCoor(dataList)
	lenDic = IORobot.obtainLength(folderName, inputName+'.fasta')

	matchThres = 8000
	nonMatchThres = 2

	separationThres = 5000
	
	newDataList= []

	for eachitem in dataList:
		name1, name2 = eachitem[-2], eachitem[-1]
		matchLen1 , matchLen2 = eachitem[4], eachitem[5]
		start1 , end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]

		if name1!= name2   and  ( min(lenDic[name1] - end1, lenDic[name2] - end2 ) > nonMatchThres or min(start1, start2) > nonMatchThres ) and matchLen1> matchThres:
			newDataList.append(eachitem)
		elif name1 == name2 and abs(start1 - start2) > separationThres  and abs(end1 - end2 ) > separationThres and matchLen1 > matchThres:
			#print eachitem
			newDataList.append(eachitem)

	#print len(newDataList), newDataList[0]
	# 351 [1, 17852, 1, 17842, 17852, 17842, 98.11, 'Segkk0', 'Segkk128']
	# assert(False)
	newDataList.sort(key = itemgetter(-2))


	if False:
		### Old method 
		repeatIntervalDic = {}
		count =0 
		for key, items in groupby(newDataList, itemgetter(-2)):
			
			listOfIntervals = []
			for eachsub in items:
				listOfIntervals.append([eachsub[0], eachsub[1]])

			if True:
				thres = 30
				B = intervalunion.intervalCover(listOfIntervals, thres)
				rangeList = intervalunion.reportMisAssemblyIntervals2(B, lenDic[key], thres, key)
				count += len(rangeList)
				if len(rangeList) > 0:
					repeatIntervalDic[key] = rangeList

		print "Count", count
		print repeatIntervalDic

	else:
		repeatIntervalDic = breakPointFinding.returnBkPtBoolSat(newDataList)

	return repeatIntervalDic

def groupCTest(folderName, repeatIntervalDic):
	'''
	Input : coveragePerContigs, repeatIntervalDic
	Output : blkDic
	'''
	json_data = open(folderName + "coveragePerContigs.json", 'r')
	coveragePerContigs = json.load(json_data)

	xDic, locDic = {}, {}
	for eachitem in repeatIntervalDic:
		xDic[eachitem], locDic[eachitem] = formatIntervalAndFilter(repeatIntervalDic[eachitem],coveragePerContigs[eachitem])
		#print xDic[eachitem], locDic[eachitem]
		
	pDic = {}
	for eachitem in xDic:
		if  len(xDic[eachitem]) > 1:
			pDic[eachitem] = groupctest.predictDecisionBoundary(xDic[eachitem], "groupCTest")
			#pDic[eachitem] = [1 for i in range(len(xDic[eachitem])-1)]

	blkDic = {}
	for eachitem in pDic:
		blkDic[eachitem] = formBlkDic(pDic[eachitem], locDic[eachitem], len(coveragePerContigs[eachitem]))

	return blkDic

def formatIntervalAndFilter2(repeatInterval, coverageData):
	tmpList = []
	initThres = 2
	lenThres = 1
	largestRange = 50000

	start, end = initThres, repeatInterval[0][0]
	end = min ( start + largestRange, end)
	tmpList.append([np.sum(coverageData[start:end]), start, end ])

	for i in range(1, len(repeatInterval)):
		start, end = repeatInterval[i-1][1], repeatInterval[i][0]
		end = min ( start + largestRange, end)
		tmpList.append([np.sum(coverageData[start:end]), start, end ])

		start, end = repeatInterval[i][0], repeatInterval[i][1]
		end = min ( start + largestRange, end)
		tmpList.append([np.sum(coverageData[start:end]), start, end ])

	start, end = repeatInterval[-1][1], len(coverageData) - initThres
	end = min ( start + largestRange, end) 
	tmpList.append([np.sum(coverageData[start:end]), start, end ])

	dataList, locList = [], []

	#print "tmpList, repeatInterval", tmpList, repeatInterval
	for eachitem in tmpList:
		if eachitem[2] - eachitem[1] > lenThres : 
			dataList.append([eachitem[0], eachitem[2] - eachitem[1]])
			locList.append([eachitem[1], eachitem[2]])

	return dataList, locList 

def formatIntervalAndFilter(repeatInterval, coverageData):
	dataList, locList = [], []
	
	initThres = 2
	lenThres = 1
	largestRange = 50000

	start, end = initThres, repeatInterval[0][0]
	if end - start > lenThres:
		locList.append([start, end])
		end = min ( start + largestRange, end)
		dataList.append([np.sum(coverageData[start:end]),  end -start ])

	for i in range(0, len(repeatInterval)-1):
		start, end = repeatInterval[i][1], repeatInterval[i+1][0]
		if end - start > lenThres:
			locList.append([start, end])
			end = min ( start + largestRange, end)
			dataList.append([np.sum(coverageData[start:end]),  end -start ])

	start, end = repeatInterval[-1][1], len(coverageData) - initThres
	if end - start > lenThres:
		locList.append([start, end])
		end = min ( start + largestRange, end)
		dataList.append([np.sum(coverageData[start:end]),  end -start ])

	return dataList, locList 

def formBlkDic(p, loc, contigLen):
	blkList = []
	blkList.append(0)

	for i in range(len(p)):
		if p[i] == 1:
			blkList.append(loc[i][1])

	blkList.append(contigLen)
	return blkList

def mainFlow(folderName):
	repeatIntervalDic = repeatFinder(folderName, "LC_filtered")
	blkDicNew = groupCTest(folderName, repeatIntervalDic)
	with open(folderName + "blkDicNew.json", 'w') as outfile:
		json.dump(blkDicNew, outfile)

#parser = argparse.ArgumentParser(description='evalmfixer')
#parser.add_argument('folderName')
#args = vars(parser.parse_args())

#mainFlow(args['folderName'])
