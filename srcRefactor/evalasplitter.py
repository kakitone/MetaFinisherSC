'''
This script is to evaluate the end to end performance in depth for A-Splitter.
Colab visualization will be accompanied.  
Done.
'''

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

print "ASplitter analysis."

def readInJSON(folderName, filename):
	json_data = open(folderName + filename, 'r')
	dataItem = json.load(json_data)
	return dataItem

def loadDataFromStages(folderName):
	stagesList = ["graphsurgery.json", "BResolution.json", "XResolution.json"]
	#stagesList = ["graphsurgery.json"]
	matchingList = []
	for filename in stagesList:
		dataList = readInJSON(folderName, filename)
		matchingList.append([dataList, filename])
	return matchingList

def isMatch(eachitem, lenDic):
	start,end = min(eachitem[2], eachitem[3]), max(eachitem[2], eachitem[3])
	isReverse = (eachitem[2] > eachitem[3])
	thres = 30
	contigName = eachitem[-1]

	if start < thres and end > lenDic[contigName] - thres and not isReverse:
		return 'f'
	elif start < thres and end > lenDic[contigName] - thres and isReverse:
		return 'r'
	else: 
		return 'n'

def findGroundTruth(folderName, mummerPath):
	# "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"

	if True:
		alignerRobot.useMummerAlign(mummerPath, folderName, "groundTruthMatch", "reference.fasta", "improved3_Double.fasta", False, "", False)
	
	dataList = alignerRobot.extractMumData(folderName, "groundTruthMatch" + "Out")

	lenDic = IORobot.obtainLength(folderName, "improved3_Double.fasta")
	lenDicRef = IORobot.obtainLength(folderName, "reference.fasta")

	dataList.sort(key=itemgetter(-1))
	
	# print "print len(lenDic), len(lenDicRef), len(dataList)", len(lenDic), len(lenDicRef), len(dataList)

	# Format of newList : [refName, refStart, refEnd, contigName]

	newList = []
	for key, items in groupby(dataList, itemgetter(-1)):
		for eachitem in items:
			if isMatch(eachitem, lenDic) == 'f':
				newList.append([eachitem[-2], eachitem[0], eachitem[1], eachitem[-1]])
				#break
			elif isMatch(eachitem, lenDic) == 'r':
				refName = eachitem[-2]
				newList.append([refName +"_r", lenDicRef[refName] - eachitem[1], lenDicRef[refName] - eachitem[1], eachitem[-1]])
				#break
	
	newList.sort()
	succList = []

	for key, items in groupby(newList, itemgetter(0)):
		tmpList = list(items)
		# print "len(tmpList)", len(tmpList)
		for i in range(len(tmpList) -1):
			succList.append([tmpList[i][-1], tmpList[i+1][-1]])

	# print "len(succList)", len(succList)
	return succList

def isInside(mySortedList, item):
	if len(mySortedList) == 0 :
		return False

	i = bisect.bisect(mySortedList,item)
	if mySortedList[max(0,i-1)] == item:
		return True
	else:
		return False

def calculatePrecisionRecall(folderName, mummerPath):
	
	matchedListAtStages = loadDataFromStages(folderName)
	allPossibleMatches = findGroundTruth(folderName,mummerPath)

	c = len(allPossibleMatches)
	alreadyUsedList = []
	allPossibleMatches.sort()
	#for eachitem in  allPossibleMatches:
	#	print eachitem
	indComponentTestList = []

	for i in range(len(matchedListAtStages)):
		alreadyUsedList.sort()
		tmpList = []
		a , b = 0, 0 
		for x in matchedListAtStages[i][0]:
			if not isInside(alreadyUsedList, x):
				tmpList.append(x)
				if isInside(allPossibleMatches, x):
					a = a +1
				else:
					b = b+1
		# print tmpList
		alreadyUsedList += tmpList

		if a>0 and c> 0:
			precision, recall = a*1.0/(a+b), a*1.0/c
		else:
			precision , recall = 1, 0
		print "i, precision, recall, TP_num, FP_num  : %d \t %f \t %f \t %d \t %d "  % (i, precision, recall, a,b)  

		indComponentTestList.append([precision, recall, a, b,  matchedListAtStages[i][1] ])

	with open(folderName + 'indComponentTestList.json', 'w') as f:
		json.dump(indComponentTestList, f)

def unitTesting():
	if False:
		kktest = [[1,2], [2,3], [3,4]]
		for item in kktest: 
			print isInside(kktest,item)
		item = [1,0]
		print isInside(kktest, item)
		item = [0,1]
		print isInside(kktest, item)
	
	if True:
		findGroundTruth("dataFolder0/")

parser = argparse.ArgumentParser(description='evalasplitter')
parser.add_argument('folderName')
parser.add_argument('mummerPath')

args = vars(parser.parse_args())

folderName = args['folderName']
mummerPath = args['mummerPath']
calculatePrecisionRecall(folderName, mummerPath)



# calculatePrecisionRecall()
# unitTesting()









