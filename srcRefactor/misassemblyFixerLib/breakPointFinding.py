import boolSat
import networkx as nx
from itertools import groupby
from operator import itemgetter

### You can even start with a trivial example and walk the way through without getting to the real code. 

### But make sure you handle the cases of reverse repeat as well 

class disjointset(object):
	def __init__(self, index):
		self.rank  = 0
		self.parent = self
		self.id = index

def union(x,y):
	xRoot = find(x)
	yRoot = find(y) 

	if xRoot == yRoot:
		return
	else:
		if xRoot.rank < yRoot.rank : 
			xRoot.parent = yRoot
		elif xRoot.rank > yRoot.rank:
			yRoot.parent = xRoot			
		else:
			yRoot.parent = xRoot
			xRoot.rank += 1

def find(x):
	if x.parent == x:
		return x
	else:
		return find(x.parent)

def formNaturalBkPts(mummerDataList):
	bkpts = []
	# dataList.append([1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])
	#  [key, repeat_index_raw, repeat_end_pt, cluster_index, contig_name, location]

	for i in range(len(mummerDataList)):
		bkpts.append([len(bkpts), i , 0, 2*i, mummerDataList[i][-2], mummerDataList[i][0] ])
		bkpts.append([len(bkpts), i , 1, 2*i + 1, mummerDataList[i][-2], mummerDataList[i][1] ])
		bkpts.append([len(bkpts), i , 2, 2*i, mummerDataList[i][-1], mummerDataList[i][2] ])
		bkpts.append([len(bkpts), i , 3, 2*i + 1, mummerDataList[i][-1], mummerDataList[i][3] ])

	return bkpts

def findRepeatDic(bkpts):
	repeatDic = {}
	for i  in range(len(bkpts)):
		eachitem = bkpts[i]
		repeatIndex, endpt  = eachitem[1], eachitem[2]
		repeatDic[str(repeatIndex) + "_" + str(endpt)] = i

	return repeatDic

def formName(i,j):
	return str(i) + "_" + str(j)

def locateEnclosedBkPts(repeatIndex, mylist, bkpts, repeatDic):
	# [key, repeat_index_raw, repeat_end_pt, cluster_index, contig_name, location]	
	# mylist [[0, 0, 0, 0, 'Contig1', 1001], [1, 0, 1, 1, 'Contig1', 2001], [2, 0, 2, 0, 'Contig2', 1001], [3, 0, 3, 1, 'Contig2', 2001]]
	# bkpts are sorted in last two coordinates 
	
	### Find Partners
	maxSize = len(bkpts) -1 

	j1  = repeatDic[formName(repeatIndex,0)]
	j2  = repeatDic[formName(repeatIndex,2)]

	# print "j1, j2", j1, j2
	j =0 	
	while j < 500:

		if j1 -j >= 0 and bkpts[j1 - j][1] == repeatIndex and bkpts[j1 - j][2] == 1:
			j1_e = j1 - j		
		if j1+j <= maxSize and bkpts[j1 + j][1] == repeatIndex and bkpts[j1 + j][2] == 1:
			j1_e = j1 + j
		j += 1 
	

	j =0 

	'''
	if  repeatIndex == 138: 
		for eachitem in bkpts[j2-10 : j2+10]:
			print eachitem
	
		print "-------"
	'''
	# print repeatIndex 
	while j < 500:
		if j2 -j >= 0 and  bkpts[j2 - j][1] == repeatIndex and bkpts[j2 - j][2] == 3:
			j2_e = j2 - j		
		if j2+j <= maxSize and bkpts[j2 + j][1] == repeatIndex and bkpts[j2 + j][2] == 3:
			j2_e = j2 + j
		j += 1 

	### Find Included Break Points
	#print "j1, j1_e, j2, j2_e", j1, j1_e, j2, j2_e
	
	involvedBkPtsList = [ bkpts[min(j1,j1_e):max(j1,j1_e)+1] , bkpts[min(j2,j2_e):max(j2,j2_e)+1] ]
	
	# print "repeatIndex", repeatIndex
	
	return involvedBkPtsList

def fillInHidden(involvedBkPtsList):
	'''
	involvedBkPtsList: 

	[
		[ 
		  [0, 0, 0, 0, 'Contig1', 1001], 
		  [1, 0, 1, 1, 'Contig1', 2001]
		], 

	    [ 
	      [2, 0, 2, 0, 'Contig2', 1001], 
	      [3, 0, 3, 1, 'Contig2', 2001]
	    ]
	]
	'''
	
	### Filter matched items
	
	firstExistDic = {}
	secondExistDic = {}
	
	involvedBkPtsList[0].sort(key = itemgetter(-2, -1))
	involvedBkPtsList[1].sort(key = itemgetter(-2, -1))
	
	for eachitem in involvedBkPtsList[0]: 
		firstExistDic[formName(eachitem[1], eachitem[3])] = True

	for eachitem in involvedBkPtsList[1]: 
		secondExistDic[formName(eachitem[1], eachitem[3])] = True
	
	newFirstList , newSecondList = [], []

	for eachitem in involvedBkPtsList[0]: 
		if not formName(eachitem[1], eachitem[3]) in secondExistDic:
			newSecondList.append(eachitem)

	for eachitem in involvedBkPtsList[1]: 
		if not formName(eachitem[1], eachitem[3]) in firstExistDic:
			newFirstList.append(eachitem)
		
	### Fill in the blanks for matched cases
	newItems = []
	# Check reverse or not 
	# Add in the missing values here 

	if involvedBkPtsList[0][0][3] == involvedBkPtsList[1][0][3]:
		reverse = False
	else:
		reverse = True

	### <--- fill in the blanks accordingly then you should be fine. 
	# [3, 0, 3, 1, 'Contig2', 2001]

	if reverse == False:
		firstBegin, secondBegin = involvedBkPtsList[0][0][-1], involvedBkPtsList[1][0][-1]

		firstContigName , secondContigName = involvedBkPtsList[0][0][-2], involvedBkPtsList[1][0][-2]

		# print "firstContigName, secondContigName", firstContigName, secondContigName

		for eachitem in newFirstList:
			newItems.append([ -1 , eachitem[1]  , -1, eachitem[3] , firstContigName ,  eachitem[-1] - secondBegin + firstBegin])

		for eachitem in newSecondList:
			newItems.append([ -1 , eachitem[1]  , -1 , eachitem[3], secondContigName, eachitem[-1] - firstBegin + secondBegin])

	else:
		firstBegin, secondBegin = involvedBkPtsList[0][0][-1], involvedBkPtsList[1][-1][-1]

		firstContigName , secondContigName = involvedBkPtsList[0][0][-2], involvedBkPtsList[1][0][-2]

		# print "reverse: firstContigName, secondContigName", firstContigName, secondContigName
		for eachitem in newFirstList:
			newItems.append([ -1 , eachitem[1]  , -1 , eachitem[3], firstContigName , secondBegin - eachitem[-1] + firstBegin])

		for eachitem in newSecondList:
			newItems.append([ -1 , eachitem[1] , -1 , eachitem[3], secondContigName, secondBegin - eachitem[-1] + firstBegin])

	'''
	print "involvedBkPtsList"
	print involvedBkPtsList[0]
	print involvedBkPtsList[1]
	print "len(newFirstList), len(newSecondList) : " , len(newFirstList), len(newSecondList)

	print len(newItems)
	'''

	return newItems

def filterDuplicate(bkptList):

	bkptList.sort(key = itemgetter(-2,-1))	
	newbkptList = []

	for key, items in groupby(bkptList, itemgetter(-2,-1)):
		newbkptList.append(list(items)[0])

	return newbkptList

def toAddPtsFormat(toAddPts, L):
	newToAddPts = []
	for i in range(len(toAddPts)):
		tmpItem = toAddPts[i]
		tmpItem[0] = L + i
		newToAddPts.append(tmpItem)

	return newToAddPts

def addHiddenBkPts(bkpts):
	# print "bkpts : ",  bkpts
	# [[0, 0, 0, 0, 'Contig1', 1001], [1, 0, 1, 1, 'Contig1', 2001], [2, 0, 2, 0, 'Contig2', 1001], [3, 0, 3, 1, 'Contig2', 2001]]
	bkPtRepeatSort = [eachitem for eachitem in bkpts]
	bkPtRepeatSort.sort(key = itemgetter(1))

	for i in range(1) :
		toAddPts = []

		bkpts.sort(key = itemgetter(-2, -1))
		repeatDic = findRepeatDic(bkpts)
		
		'''
		print "bkpts"
		
		for i in range(len(bkpts)):
			print i , bkpts[i]
		'''

		for key , items in groupby(bkPtRepeatSort, itemgetter(1)):
			repeatIndex = key 
			mylist = list(items)
			involvedBkPts = locateEnclosedBkPts(repeatIndex, mylist, bkpts, repeatDic)
			newItems = fillInHidden(involvedBkPts)
			toAddPts += newItems
		
		'''
		toAddPts = toAddPtsFormat(toAddPts, len(bkpts))
		bkpts = toAddPts + bkpts
		'''

		toAddPts = toAddPtsFormat(toAddPts, len(bkpts))
		newBkpts =  toAddPts + bkpts

		bkptsfiltered = filterDuplicate(newBkpts)

	return bkptsfiltered

def clusterBkPts(newbkts):	
	### Get merge pairs 
	#newbkts = oldbkts
	#newbkts = filterDuplicate(oldbkts)
	newbkts.sort(key = itemgetter(-2, -1))
	mergingPair = []
	basicItems = []

	for key , items in groupby(newbkts, itemgetter(-2)):
		mylist = list(items)
		mylist.sort(key = itemgetter(-1))

		basicItems.append(mylist[0][3])
		for i in range(len(mylist)-1):
			basicItems.append(mylist[i+1][3])

			if abs(mylist[i+1][-1] - mylist[i][-1]) < 60:
				mergingPair.append([mylist[i+1][3], mylist[i][3]])

	### Get Use disjoint set union
	distinctBasic = []
	for eachitem in set(basicItems):
		distinctBasic.append(eachitem)

	totalClusters = [disjointset(i) for i in distinctBasic]
	key2IndexMap = { distinctBasic[i] : i for i in range(len(distinctBasic)) }

	#print key2IndexMap
	#print len(totalClusters)

	for eachpair in mergingPair:
		x, y = totalClusters[key2IndexMap[eachpair[0]]], totalClusters[key2IndexMap[eachpair[1]]]
		union(x, y)

	### Re write the newbks to form newClusteredBks
	### Something is wrong here 

	newClusteredBks = [] 
	newbkts.sort(key = itemgetter(-2, -1))

	for eachitem in newbkts:		
		tmpitem = eachitem
		tmpitem[3] = find(totalClusters[key2IndexMap[tmpitem[3]]]).id
		newClusteredBks.append(tmpitem)

	newClusteredBks.sort(key = itemgetter(-2, -1))
	newClusteredBks2 = [ newClusteredBks[0] ]

	for i in range(1, len(newClusteredBks)):
		if  newClusteredBks[i][3] == newClusteredBks[i-1][3] and newClusteredBks[i][-2] == newClusteredBks[i-1][-2] and abs(newClusteredBks[i][-1] - newClusteredBks[i-1][-1] ) < 60:
			assert(True)
		else: 
			newClusteredBks2.append(newClusteredBks[i])

	return newClusteredBks2

def findSeqs(newClusteredBks):
	linearSeqs, bkPt2ClustersMapDic = [], {}
	newClusteredBks.sort(key = itemgetter(-2, -1))
	for key, items in groupby(newClusteredBks, itemgetter(-2)):
		tmpSeq = []
		for eachitem in items:
			tmpSeq.append(eachitem[0])
			bkPt2ClustersMapDic[eachitem[0]] = eachitem[3]
		linearSeqs.append(tmpSeq)

	return linearSeqs, bkPt2ClustersMapDic

def formClusterName(myList, bkPt2ClustersMapDic):
	newList = []
	for eachitem in myList:
		newList.append(bkPt2ClustersMapDic[eachitem])
	return newList 

def filterList(myList):
	mynewlist = []
	myList.sort()
	for key, items in groupby(myList):
		mynewlist.append(key)

	return mynewlist

def findInOutList(linearSeqs, bkPt2ClustersMapDic):
	'''
	This may be the easiest to beging
	
	Input : 	linearSeqs  =   [[1 , 2, 3, 4, 5, 6], [7,8,9,10], [11, 12, 13] ] <--- cannot duplicated , bkPt2ClustersMapDic   
	Output : 	setOfChoices = [[[11], [12] ],  [ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] ]
	Algorithm : 
		1. Transform to clusterForm 		
		2. Scan for Common pairs 
		3. Find Neighboring items 
		4. Report a pair then 

	'''
	endAppendedList = []
	clusterFilledList = []
	indicesPairList = []
	setOfChoices = []

	### Forming endAppendedList
	maxIndex = -1 
	for eachSeq in linearSeqs:
		for eachitem in eachSeq :
			if eachitem > maxIndex : 
				maxIndex = eachitem 

	maxIndex += 1


	for eachSeq in linearSeqs : 
		
		bkPt2ClustersMapDic[maxIndex] = maxIndex
		tmpSeq = [maxIndex ]	
		maxIndex += 1

		for eachitem in eachSeq :
			tmpSeq.append(eachitem)

		bkPt2ClustersMapDic[maxIndex] = maxIndex
		tmpSeq.append(maxIndex)
		maxIndex += 1

		endAppendedList.append(tmpSeq)

	### Forming clusterFilledList

	for eachSeq in endAppendedList:
		tmpSeq = []
		for eachitem in eachSeq : 
			tmpSeq.append(bkPt2ClustersMapDic[eachitem])
		clusterFilledList.append(tmpSeq)


	### Forming indicesPairList
	### Keep track of the end points 
	for i in range(len(clusterFilledList)):
		for j in range(1, len(clusterFilledList[i]) -1):
			indicesPairList.append([clusterFilledList[i][j-1], clusterFilledList[i][j], i, j])
	
	indicesPairList.sort(key = itemgetter(0,1))	

	### Forming indicesPairList
	for key, items in groupby(indicesPairList, itemgetter(0,1)):
		inList, outList, prevList, nextList = [], [] , [], []

		for eachitem in items:
			i, j = eachitem[-2], eachitem[-1]
			inItem = endAppendedList[i][j-1]
			outItem = endAppendedList[i][j]

			prevList.append(clusterFilledList[i][j-2])
			nextList.append(clusterFilledList[i][j+1])

			inList.append(inItem)
			outList.append(outItem)

		inList = filterList(inList)
		outList = filterList(outList)

		prevList = filterList(prevList)
		nextList = filterList(nextList)

		if len(prevList) > 1 and len(nextList) > 1 :
			setOfChoices.append([inList, outList])


	return setOfChoices 

def returnBkPts(combineCutList, newClusteredBks):
	bkPtDic = {}
	
	for eachitem in  newClusteredBks:
		bkPtDic[eachitem[0]] = [eachitem[-2], eachitem[-1]]

	returnFormattedBkPtList = []

	for eachitem in combineCutList:
		returnFormattedBkPtList.append(bkPtDic[eachitem]) 

	return returnFormattedBkPtList

def pipelineOriginalMethod(returnFormattedBkPtList):
	### May change a bit format depending on need. 
	returnFormattedBkPtList.sort()
	repeatIntervalDic = {}
	# {'Segkk1': [[2500001, 2500016], [2511985, 2512000]], 'Segkk0': [[2500001, 2500016], [2511985, 2512000]]}

	for key, items in groupby(returnFormattedBkPtList, itemgetter(0)):
		repeatIntervalDic[key] = []
		for eachitem in items:
			repeatIntervalDic[key].append([ eachitem[-1] , eachitem[-1] + 15 ] )

	return repeatIntervalDic

def returnBkPtBoolSat(mummerDataList):

	bkpts = formNaturalBkPts(mummerDataList)
	print "len(bkpts)", len(bkpts)

	newbkts = addHiddenBkPts(bkpts)
	print "len(newbkts)", len(newbkts)

	newClusteredBks = clusterBkPts(newbkts)
	print "len(newClusteredBks)", len(newClusteredBks)

	linearSeqs, bkPt2ClustersMapDic = findSeqs(newClusteredBks)
	print "len(linearSeqs)", len(linearSeqs)
	print "len(bkPt2ClustersMapDic)", len(bkPt2ClustersMapDic)

	setOfChoices = findInOutList(linearSeqs, bkPt2ClustersMapDic)
	print "len(setOfChoices)", len(setOfChoices)

	combineCutList = boolSat.findMinSat(setOfChoices)
	print "len(combineCutList)", len(combineCutList)

	returnFormattedBkPtList = returnBkPts(combineCutList, newClusteredBks)
	print "len(returnFormattedBkPtList)", len(returnFormattedBkPtList)

	repeatIntervalDic = pipelineOriginalMethod(returnFormattedBkPtList)
	print "len(repeatIntervalDic)", len(repeatIntervalDic)

	#assert(False)
	return repeatIntervalDic







