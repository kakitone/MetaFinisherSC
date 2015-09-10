from ..repeatPhaserLib.finisherSCCoreLib import IORobot
from ..repeatPhaserLib.finisherSCCoreLib import alignerRobot
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper 
import bisect
import adaptorFix
import groupctest
import intervalunion

from operator import itemgetter
from itertools import groupby
import json
import numpy as np


def repeatFinder(folderName, inputName):
	'''
	Input : dataList
	Output : repeatIntervalDic
	'''
	dataList = alignerRobot.extractMumData(folderName, "self"+inputName+"Out")
	dataList = alignerRobot.transformCoor(dataList)
	lenDic = IORobot.obtainLength(folderName, inputName+'.fasta')

	matchThres = 10000
	nonMatchThres = 500
	
	newDataList= []

	for eachitem in dataList:
		name1, name2 = eachitem[-2], eachitem[-1]
		matchLen1 , matchLen2 = eachitem[4], eachitem[5]
		start1 , end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]

		if name1!= name2   and  ( min(lenDic[name1] - end1, lenDic[name2] - end2 ) > nonMatchThres or min(start1, start2) > nonMatchThres ) and matchLen1> matchThres:
			newDataList.append(eachitem)

	print len(newDataList), newDataList[0]
	# 351 [1, 17852, 1, 17842, 17852, 17842, 98.11, 'Segkk0', 'Segkk128']
	newDataList.sort(key = itemgetter(-2))

	repeatIntervalDic = {}

	for key, items in groupby(newDataList, itemgetter(-2)):
		
		listOfIntervals = []
		for eachsub in items:
			listOfIntervals.append([eachsub[0], eachsub[1]])

		repeatIntervalDic[key] = intervalunion.intervalUnion(listOfIntervals) 

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
			#print xDic[eachitem]
			pDic[eachitem] = groupctest.predictDecisionBoundary(xDic[eachitem], "groupCTest")
			print xDic[eachitem], pDic[eachitem]
			
	blkDic = {}
	for eachitem in pDic:
		blkDic[eachitem] = formBlkDic(pDic[eachitem], locDic[eachitem], len(coveragePerContigs[eachitem]))

	return blkDic


def formatIntervalAndFilter(repeatInterval, coverageData):
	tmpList = []
	initThres = 500
	lenThres = 1000

	start, end = initThres, repeatInterval[0][0]
	# Potential hacks cut head/end too

	tmpList.append([np.sum(coverageData[start:end]), start, end ])

	for i in range(0, len(repeatInterval)-1):
		start, end = repeatInterval[i][1], repeatInterval[i+1][0]
		tmpList.append([np.sum(coverageData[start:end]), start, end ])

	start, end = repeatInterval[-1][1], len(coverageData) - initThres
	tmpList.append([np.sum(coverageData[start:end]), start, end ])

	dataList, locList = [], []

	#print "tmpList, repeatInterval", tmpList, repeatInterval
	for eachitem in tmpList:
		if eachitem[2] - eachitem[1] > lenThres : 
			dataList.append([eachitem[0], eachitem[2] - eachitem[1]])
			locList.append([eachitem[1], eachitem[2]])

	return dataList, locList 


def formBlkDic(p, loc, contigLen):
	blkList = []
	blkList.append(0)
	for i in range(len(p)):
		if p[i] == 1:
			blkList.append(loc[i][0])

	blkList.append(contigLen)
	return blkList


def unitTest():
	folderName = "dataFolder/"
	repeatIntervalDic = repeatFinder(folderName, "LC_filtered")
	if False:
		for eachitem in repeatIntervalDic:
			print repeatIntervalDic[eachitem]

	blkDicNew = groupCTest(folderName, repeatIntervalDic)

	with open(folderName + "blkDicNew.json", 'w') as outfile:
		json.dump(blkDicNew, outfile)

unitTest()










