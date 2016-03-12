'''
Input : [ [ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] , [ ... ] , [] ]
Output : [ 4, 5, 6, 9, 10, ...]
Algorithm :
	1. Find Related clusters
	2. For each clusters : 
		if len(clusters) <= K 
			Exhaust
		else
			Greedily choose/randomly pick

'''
import networkx as nx
import itertools

inputList = []
ouputList = []

def overlapTrue(myList1, myList2):
	combine1 = []
	combine2 = []
	for eachitem in myList1:
		combine1 += eachitem 

	for eachitem in myList2:
		combine2 += eachitem 
	
	if len(set(combine1).intersection(set(combine2))) > 0:
		return True
	else:
		return False

def findRelatedClusters(inputList):
	clusters = []
	G=nx.Graph()
	for i in range(len(inputList)):
		G.add_node(i)

	for i in range(len(inputList)):
		for j in range(len(inputList)):
			if i != j :
				if overlapTrue(inputList[i] , inputList[j]):
					G.add_edge(i,j)

	G.to_undirected()

	clusters =  [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]

	return clusters

def totalScore(listTuple):
	score = 0
	myList = []
	for eachitem in listTuple:
		myList += eachitem 

	score = len(set(myList))
	return score

def exhaustList(myList):
	'''
	Permute all possible combinations 
	Output a score and report the highest one 
	'''
	
	minScore = 10**5
	minElem = []
	for element in itertools.product(*myList):
		myScore =  totalScore(element) 
		if myScore <= minScore:
			minElem = element
			minScore = myScore

	minList = [] 

	for eachitem in minElem:
		minList += eachitem
	
	return [ eachitem for eachitem in set(minList)]

def greedyPickList(myList):
	newList = []
	for eachitem in myList:
		newList += eachitem[0]
	return [ eachitem for eachitem in set(newList)]

def formSubList(inputList, eachcluster):
	return [inputList[i] for i in eachcluster]

def findMinSat(inputList):
	K = 5 
	clusters = findRelatedClusters(inputList)
	combineCutList = []

	for eachcluster in clusters:
		subList = formSubList(inputList, eachcluster)
		
		if len(eachcluster) <= K:
			combineCutList += exhaustList(subList)
		else:
			combineCutList += greedyPickList(subList)

	return combineCutList

def testing():
	if False:
		inputList =  [[[11], [12] ],  [ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] ] 
		print findRelatedClusters(inputList)

	if False:
		exhaustList([[ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] ])

	inputList = [[[11], [12] ],  [ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] ]
	print findMinSat(inputList)

# testing()
# findMinSat(inputList)

