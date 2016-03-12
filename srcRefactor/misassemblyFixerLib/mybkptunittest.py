import os
import unittest
import breakPointFinding
import boolSat
from operator import itemgetter

class bkFindingAlgoTest(unittest.TestCase):    

	def setUp(self):
		print "Set up : Started"

	def test_Test1(self):
		dataList = []
		dataList.append([1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])

		expectedOutput = [ ['Contig1', 1001] , ['Contig2', 1001] ]

	def test_Test2(self):
		dataList = []
		dataList.append([1001, 2001, 2001, 1001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])

		expectedOutput = [ ['Contig1', 1001] , ['Contig2', 2001] ]

	def test_Test3(self):
		dataList = []
		dataList.append([1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])
		dataList.append([1501, 2501, 1501, 2501, 1000, 1000, 100.0, 8000, 8000, 'Contig2', 'Contig3'])

		expectedOutput = [ ['Contig1', 1501] , ['Contig2', 1501], ['Contig3', 1501] ]

	def test_Test4(self):
		dataList = []
		dataList.append([1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])
		dataList.append([1501, 2501, 2501, 1501, 1000, 1000, 100.0, 8000, 8000, 'Contig2', 'Contig3'])

		expectedOutput = [ ['Contig1', 1501] , ['Contig2', 1501], ['Contig3', 2501] ]

	def test_findInOutList(self):
		'''	
		Input : 	linearSeqs  =   [[1 , 2, 3, 4, 5, 6], [7,8,9,10], [11, 12, 13] ] <--- cannot duplicated , bkPt2ClustersMapDic   
		Output : 	setOfChoices = [[[11], [12] ],  [ [1,2,3] , [4,5,6] ] , [ [ 4, 5], [ 9, 10] ] ]
		'''

		linearSeqs, bkPt2ClustersMapDic = [ [1, 2 ], [3, 4] ], {1:1, 2:2, 3:1, 4:2  }
		setOfChoices = breakPointFinding.findInOutList(linearSeqs, bkPt2ClustersMapDic)
		setOfChoices.sort() 				


		expectedSetOfChoices = [[[1, 3], [2, 4]]]
 		
		assert(setOfChoices == expectedSetOfChoices)	

		linearSeqs, bkPt2ClustersMapDic = [ [1, 2 , 3 ], [4, 5, 6, 7], [8,9,10] ], { 1:1 ,2:2 , 3:3 , 4: 1, 5:2, 6:3, 7:4, 8:2, 9:3, 10:4 }
		setOfChoices = breakPointFinding.findInOutList(linearSeqs, bkPt2ClustersMapDic)
		setOfChoices.sort()

		expectedSetOfChoices = [[[2, 5, 8], [3, 6, 9]]]
		
		assert(setOfChoices == expectedSetOfChoices)

	def test_formNaturalBkPts(self):
		mummerDataList = []
		mummerDataList.append([1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])

		bkpts = breakPointFinding.formNaturalBkPts(mummerDataList)

		expectedBkPtsList = []
		### [key, repeat_index_raw, repeat_end_pt, cluster_index, contig_name, location]
		
		expectedBkPtsList.append([0, 0, 0, 0, "Contig1", 1001])
		expectedBkPtsList.append([1, 0, 1, 1, "Contig1", 2001])
		expectedBkPtsList.append([2, 0, 2, 0, "Contig2", 1001])
		expectedBkPtsList.append([3, 0, 3, 1, "Contig2", 2001])

		assert(bkpts == expectedBkPtsList)

		mummerDataList = []
		mummerDataList.append([1001, 2001, 2001, 1001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2'])
		bkpts = breakPointFinding.formNaturalBkPts(mummerDataList)
		
		expectedBkPtsList = [[0, 0, 0, 0, 'Contig1', 1001], \
							 [1, 0, 1, 1, 'Contig1', 2001], \
							 [2, 0, 2, 0, 'Contig2', 2001], \
							 [3, 0, 3, 1, 'Contig2', 1001]]

 		assert(bkpts == expectedBkPtsList)

	def test_addHiddenBkPts(self):
		bkpts = []
				
		bkpts.append([0, 0, 0, 0, "Contig1", 1001])
		bkpts.append([1, 0, 1, 1, "Contig1", 2001])
		bkpts.append([2, 0, 2, 0, "Contig2", 1001])
		bkpts.append([3, 0, 3, 1, "Contig2", 2001])

		newbkts = breakPointFinding.addHiddenBkPts(bkpts)		
		
		expectedNewbkts = [[0, 0, 0, 0, 'Contig1', 1001], \
						   [1, 0, 1, 1, 'Contig1', 2001], \
						   [2, 0, 2, 0, 'Contig2', 1001], \
						   [3, 0, 3, 1, 'Contig2', 2001]]

		assert( newbkts == expectedNewbkts)

		# print "--------------------------------"
		'''
		[1001, 2001, 1001, 2001, 1000, 1000, 100.0, 8000, 8000, 'Contig1', 'Contig2']

		[1501, 2501, 2501, 1501, 1000, 1000, 100.0, 8000, 8000, 'Contig2', 'Contig3']
		'''
		bkpts = []
				
		bkpts.append([0, 0, 0, 0, "Contig1", 1001])
		bkpts.append([1, 0, 1, 1, "Contig1", 2001])
		bkpts.append([2, 0, 2, 0, "Contig2", 1001])
		bkpts.append([3, 0, 3, 1, "Contig2", 2001])
		
		bkpts.append([4, 1, 0, 2, "Contig2", 1501])
		bkpts.append([5, 1, 1, 3, "Contig2", 2501])
		bkpts.append([6, 1, 3, 2, "Contig3", 2501])
		bkpts.append([7, 1, 2, 3, "Contig3", 1501])


		newbkts = breakPointFinding.addHiddenBkPts(bkpts)		
		newbkts.sort(key = itemgetter(-2, -1))
		
		#for eachitem in  newbkts:
		#	print eachitem 

	def test_clusterBkPts(self):
		newbkts  = [	\
						[0, 0, 0, 0, 'Contig1', 1001], \
						[8, 1, -1, 2, 'Contig1', 1501],\
						[1, 0, 1, 1, 'Contig1', 2001],\
						[2, 0, 2, 0, 'Contig2', 1001],\
						[4, 1, 0, 2, 'Contig2', 1501],\
						[3, 0, 3, 1, 'Contig2', 2001],\
						[5, 1, 1, 3, 'Contig2', 2501],\
						[7, 1, 2, 3, 'Contig3', 1501],\
						[9, 0, -1, 1, 'Contig3', 2001],\
						[10, 0, -1, 4, 'Contig3', 2002],\
						[6, 1, 3, 2, 'Contig3', 2501]\
				   ]

		newClusteredBks = breakPointFinding.clusterBkPts(newbkts)

		expectedNewClusteredBks = [ \
			[0, 0, 0, 0, 'Contig1', 1001], \
			[8, 1, -1, 2, 'Contig1', 1501],\
			[1, 0, 1, 4, 'Contig1', 2001],\
			[2, 0, 2, 0, 'Contig2', 1001],\
			[4, 1, 0, 2, 'Contig2', 1501],\
			[3, 0, 3, 4, 'Contig2', 2001],\
			[5, 1, 1, 3, 'Contig2', 2501],\
			[7, 1, 2, 3, 'Contig3', 1501],\
			[9, 0, -1, 4, 'Contig3', 2001],\
			[6, 1, 3, 2, 'Contig3', 2501] \
		] 

		assert(newClusteredBks == expectedNewClusteredBks)

	def test_findSeqs(self):
		newClusteredBks = [ \
			[0, 0, 0, 0, 'Contig1', 1001], \
			[8, 1, -1, 2, 'Contig1', 1501],\
			[1, 0, 1, 4, 'Contig1', 2001],\
			[2, 0, 2, 0, 'Contig2', 1001],\
			[4, 1, 0, 2, 'Contig2', 1501],\
			[3, 0, 3, 4, 'Contig2', 2001],\
			[5, 1, 1, 3, 'Contig2', 2501],\
			[7, 1, 2, 3, 'Contig3', 1501],\
			[9, 0, -1, 4, 'Contig3', 2001],\
			[6, 1, 3, 2, 'Contig3', 2501] \
		] 

		linearSeqs, bkPt2ClustersMapDic = breakPointFinding.findSeqs(newClusteredBks)

		expectedLinearSeqs = [[0, 8, 1], [2, 4, 3, 5], [7, 9, 6]]
		expectedBkPt2ClustersMapDic = {0: 0, 1: 4, 2: 0, 3: 4, 4: 2, 5: 3, 6: 2, 7: 3, 8: 2, 9: 4}

		assert(linearSeqs == expectedLinearSeqs)
		assert(bkPt2ClustersMapDic == expectedBkPt2ClustersMapDic )

	def test_boolSat(self):

		inputList = [ \
					  [ [11], [12] ], \
					  [ [1,2,3] , [4,5,6] ] , \
					  [ [ 4, 5], [ 9, 10] ] \
					]

		combineCutList = boolSat.findMinSat(inputList)
		expectedCombineCutList = [4, 5, 6, 12]

		assert(combineCutList == expectedCombineCutList)

	def test_returnBkPts(self):
		'''
		Estimate time : 2 hours 

		1. Need to hack into the existing system to see what is being expected as an input here and then generate accordingly. 
		
		2. Afterwards, it is time to submit a job for this one and see how things are going 
		
		3. May need a test for the boolSat. 
		
		'''
		combineCutList = [8, 4, 6]
		
		newClusteredBks = [ \
			[0, 0, 0, 0,  'Contig1', 1001], \
			[8, 1, -1, 2, 'Contig1', 1501],\
			[1, 0, 1, 4, 'Contig1', 2001], \
			[2, 0, 2, 0, 'Contig2', 1001], \
			[4, 1, 0, 2, 'Contig2', 1501], \
			[3, 0, 3, 4, 'Contig2', 2001], \
			[5, 1, 1, 3, 'Contig2', 2501], \
			[7, 1, 2, 3, 'Contig3', 1501], \
			[9, 0, -1, 4, 'Contig3', 2001], \
			[6, 1, 3, 2, 'Contig3', 2501] \
		] 

		returnFormattedBkPtList = breakPointFinding.returnBkPts(combineCutList, newClusteredBks)
		
		print returnFormattedBkPtList

		assert(True)

	def tearDown(self):
		print "Teardown : Started"

def main():
    unittest.main()
    
if __name__ == '__main__':
    main()










