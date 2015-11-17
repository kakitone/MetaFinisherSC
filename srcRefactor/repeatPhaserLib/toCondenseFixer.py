from itertools import groupby
from operator import itemgetter

import abunHouseKeeper
import abunGraphLib


from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import alignerRobot

import json


thresMiddleContig = abunHouseKeeper.abunGlobalSplitParameterRobot.thresMiddleContig
lenthresMiddleContig = abunHouseKeeper.abunGlobalSplitParameterRobot.lenthresMiddleContig
bestMatchContigOnly = abunHouseKeeper.abunGlobalSplitParameterRobot.bestMatchContigOnly


# Make sure G is of type seqGraphDynamic so as to use all the library
def noGoZoneDefiner(G, folderName):
	print "noGoZoneDefiner"
	noPrevList, noNextList = findNoHeads(G)
	print "len(noPrevList), len(noNextList)", len(noPrevList), len(noNextList)
	

	noGoPrev, noGoNext = findNoGoZone(noPrevList, noNextList, folderName)
	print "len(noGoPrev), len(noGoNext)", len(noGoPrev), len(noGoNext)
	
	for eachitem in noGoPrev:
		print "noGoPrev", eachitem

	for eachitem in noGoNext:
		print "noGoNext", eachitem
	

	with open( folderName + "noGoPrev.json", 'w') as f:
		json.dump(noGoPrev, f)    


	with open( folderName + "noGoNext.json", 'w') as f:
		json.dump(noGoNext, f)    

	#assert(1==2)


	#Gnew = removeEdgeInBatch(G, noGoPrev, noGoNext)

	return G

def findNoHeads(G):
	G.findStartEndList()
	return 	G.myStartList, G.myEndList 

def findNoGoZone(noPrevList, noNextList, folderName):
	noGoPrev, noGoNext = [] , []
	
	noGoPrev = findNoGoByNoHeads(noPrevList, 'L', folderName) 
	
	noGoNext = findNoGoByNoHeads(noNextList, 'R', folderName)

	return noGoPrev, noGoNext

def findNoGoByNoHeads(noGoList, side, folderName):
	noGoListNew = []

	sortedContigList,  sortedReadList, sortedContigDic, sortedReadDic =\
		formSortedDataList(folderName)


	lenDicContig = IORobot.obtainLength(folderName, "mFixed_Double.fasta" )
	lenDicRead = IORobot.obtainLength(folderName, "phasingSeedName_Double.fasta")

	for x in noGoList:
		rList = findAttachedReads(x, side, folderName,sortedContigList,sortedContigDic, lenDicContig,lenDicRead)
		cList = findAttachedContigs(rList, side, folderName, sortedReadList, sortedReadDic, lenDicContig,lenDicRead)

		if bestMatchContigOnly == False:
			bestContigIDList = findBreakContigAdv(cList)
		else:
			bestContigIDList = findBreakContig(cList)

		if len(rList) > 0 and len(cList) > 0:
			print "x, side, len(rList), len(cList), len(bestContigIDList)",\
				 abunHouseKeeper.parseIDToName(x,'C',0), side, len(rList), len(cList), len(bestContigIDList)
			print "cList", bestContigIDList
 
		noGoListNew = noGoListNew + bestContigIDList


	return noGoListNew

def formSortedDataList(folderName):
	sortedContigList,  sortedReadList, sortedContigDic, sortedReadDic =\
		[], [] , {}, {}

	dataList = alignerRobot.extractMumData(folderName, "phaseStringCROut")

	sortedContigList = sorted(dataList, key = itemgetter(-2))

	sortedContigDic[sortedContigList[0][-2]] = 0
	for i in range(1, len(sortedContigList)):
		if sortedContigList[i][-2] != sortedContigList[i-1][-2]:
			sortedContigDic[sortedContigList[i][-2]] = i



	sortedReadList = sorted(dataList, key = itemgetter(-1))
	sortedReadDic[sortedReadList[0][-1]] = 0

	for i in range(1, len(sortedReadList)):
		if sortedReadList[i][-1] != sortedReadList[i-1][-1]:
			sortedReadDic[sortedReadList[i][-1]] = i

	return sortedContigList,  sortedReadList, sortedContigDic, sortedReadDic

def removeEdgeInBatch(G, noGoPrev, noGoNext):

	for v in noGoPrev:
		G.clearIn(abunHouseKeeper.parseEdgeNameToID(v, 'C'))

	for v in noGoNext:
		G.clearOut(abunHouseKeeper.parseEdgeNameToID(v, 'C'))
	
	return G 

def findAttachedReads(x, side, folderName,sortedContigList,sortedContigDic, lenDicContig,lenDicRead):
	rList = [] 

	
	'''
	Format : 
	  [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [ IDY]  | [TAGS]
	=====================================================================================
       1      562  |      819     1418  |      562      600  |    84.72  | Contig0_d	Read121_d
       1      562  |     4077     3478  |      562      600  |    84.72  | Contig0_d	Read121_p
       1      564  |      656       68  |      564      589  |    90.13  | Contig0_d	Read382_d
       1      564  |     6996     7584  |      564      589  |    90.13  | Contig0_d	Read382_p
       1      571  |     1386      815  |      571      572  |    86.60  | Contig0_d	Read421_d

	'''


	thres = thresMiddleContig

	key = abunHouseKeeper.parseIDToName(x, 'C', 0)
	if key in sortedContigDic:
		tmp = sortedContigDic[key]
		
		while tmp < len(sortedContigList) and sortedContigList[tmp][-2] == key:
			eachsub = sortedContigList[tmp]
			if overlapCR(eachsub, side, thres, lenDicContig,lenDicRead):
				rList.append(eachsub[-1])

			tmp = tmp + 1


		distinctRList = abunHouseKeeper.getDistinct(rList)

	else:
		distinctRList = []
	return distinctRList

def overlapCR(eachsub, side, thres, lenDicContig,lenDicRead):
	if side == 'R':
		return THMatch([eachsub[0], eachsub[1], lenDicContig[eachsub[-2]]], [eachsub[2], eachsub[3], lenDicRead[eachsub[-1]]], thres, min(eachsub[4:6]) )
	else:
		return THMatch( [eachsub[2], eachsub[3], lenDicRead[eachsub[-1]]],[eachsub[0], eachsub[1], lenDicContig[eachsub[-2]]],thres, min(eachsub[4:6])   )

def overlapCRJustREnd(eachsub, side, thres, lenDicContig,lenDicRead):
	if side == 'R':
		return THMatchJustREnd([eachsub[0], eachsub[1], lenDicContig[eachsub[-2]]], [eachsub[2], eachsub[3], lenDicRead[eachsub[-1]]], thres, side, min(eachsub[4:6])  )
	else:
		return THMatchJustREnd( [eachsub[2], eachsub[3], lenDicRead[eachsub[-1]]],[eachsub[0], eachsub[1], lenDicContig[eachsub[-2]]],thres , side , min(eachsub[4:6]) )

def THMatchJustREnd(seg1Info, seg2Info, thres , side, matchLen):
	seg1Start, seg1End, seg1Len = seg1Info[0] ,seg1Info[1],seg1Info[2]
	seg2Start, seg2End, seg2Len = seg2Info[0] ,seg2Info[1],seg2Info[2]

	if side == 'R':
		if  seg2End > seg2Len - thres and matchLen > lenthresMiddleContig:
			return True	
		else:
			return False
	else : 
		if seg1Start < thres and matchLen > lenthresMiddleContig:
			return True	
		else:
			return False

def THMatch(seg1Info, seg2Info, thres, matchLen ):
	seg1Start, seg1End, seg1Len = seg1Info[0] ,seg1Info[1],seg1Info[2]
	seg2Start, seg2End, seg2Len = seg2Info[0] ,seg2Info[1],seg2Info[2]

	if seg1End > seg1Len - thres and seg2Start < thres  and matchLen > lenthresMiddleContig:
		return True	
	else:
		return False

def findAttachedContigs(rList, side, folderName, sortedReadList, sortedReadDic, lenDicContig,lenDicRead):
	cList = [] 

	thres = thresMiddleContig
	for r in rList:
		if r in sortedReadDic:
			tmp = sortedReadDic[r]
			while tmp < len(sortedReadList) and sortedReadList[tmp][-1] == r:	
				eachsub = sortedReadList[tmp]
				if overlapCRJustREnd(eachsub, side, thres,  lenDicContig, lenDicRead):
					cList.append([eachsub[-2], r])
				tmp = tmp + 1

	newCList = abunHouseKeeper.getDistinct(cList)

	return newCList

def findBreakContig(cList):

	#print "len(cList)" , len(cList)
	bestContigID = [] 
	cList.sort()
	newCList = []


	for key, items in groupby(cList, itemgetter(0)):
		ct = len(list(items))
		newCList.append([ct, key])
	newCList.sort(reverse= True)

	if len(newCList) > 0:
		bestContigID = newCList[0][1]

		return 	[bestContigID]
	else:
		return []


def findBreakContigAdv(cList):

	#print "len(cList)" , len(cList)
	bestContigID = [] 
	cList.sort()
	keyList = []
	
	for key, items in groupby(cList, itemgetter(0)):
		ct = len(list(items))
		keyList.append(key)

	if len(keyList) > 0:
		return keyList	
	else:
		return []

