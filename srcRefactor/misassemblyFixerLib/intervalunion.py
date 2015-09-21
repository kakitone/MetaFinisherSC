def intervalUnion(xList):
	B = []
	xList.sort()

	tmpL , tmpR = xList[0][0], xList[0][1]

	for i in range(1, len(xList)):
		curL, curR = xList[i][0], xList[i][1]
		if curL < tmpR:
			if curR <= tmpR:
				pass
			else :
				tmpR = curR
		else:
			B.append([tmpL, tmpR])
			tmpL, tmpR = curL, curR 

	B.append([tmpL, tmpR])

	return B

def intervalCover(xList, thres):
	### Find minimal exact cover 	
	B = findMinimalExactCover(xList)
	#print "B", B
	### Find minimal approximate cover 
	B2 =  findMinimalApproxCover(B, thres)
	#print "B2", B2

	return B2

def numeric_compare(x, y):
	if x[0] - y[0] != 0:
		return 1 if x[0] - y[0] > 0 else -1
	else:
		return 1 if y[1] - x[1] > 0 else -1

def findMinimalExactCover(xList):
	xList.sort(cmp = numeric_compare)
	Bnew = [] 
	Bnew = [xList[0]]
	ck = False

	tmpL, tmpR = xList[0][0], xList[0][1]
	tmpExtend = [tmpL, tmpR]

	for i in range(1, len(xList)):
		curL, curR = xList[i][0], xList[i][1]

		if curL < tmpR:
			if curR > tmpExtend[1]:
				tmpExtend = [curL, curR]
				ck = True
		else:
			if ck:
				Bnew.append(tmpExtend)
				ck = False
			
			Bnew.append([curL, curR])
			tmpL, tmpR= curL, curR
			tmpExtend = [curL, curR]

	if ck :
		Bnew.append(tmpExtend)

	return Bnew

def findMinimalApproxCover(xList, t = 30):
	Bnew = []
	prevInterval =[0, 0]

	for i in range(len(xList)):
		if i == len(xList) - 1:
			if isInsideInterval(xList[i],prevInterval, prevInterval, t):
				pass
			else:
				Bnew.append(xList[i])
		else:
			if isInsideInterval(xList[i], prevInterval, xList[i+1], t):
				pass
			else:
				Bnew.append(xList[i])
				prevInterval = xList[i]
	return  Bnew

def isInsideInterval(curInterval,A, B, t):
	check = False
	#print A, B, t
	if A[0] - t < curInterval[0] < curInterval[1] < A[1] + t:
		check = True
	if B[0] - t < curInterval[0] < curInterval[1] < B[1] + t:
		check = True
	if A[0] - t < curInterval[0] < curInterval[1] < B[1] + t and A[1] + t > B[0] - t:
		check = True

	return check 

def reportMisAssemblyIntervals(xList, lencontig, t = 30):
	mList = []
	
	if xList[0][0] > t: 
		mList.append([xList[0][0]-t/2, xList[0][0]+t/2]) 

	for i in range(len(xList) -1):
		diff = xList[i+1][0] - xList[i][1]
		if  diff < 0:
			mList.append([xList[i+1][0],xList[i][1] ])
		elif diff > t:
			mList.append([xList[i][1] - t/2, xList[i][1] + t/2])
			mList.append([xList[i+1][0] - t/2, xList[i+1][0] + t/2])
		else : 
			mList.append([xList[i][1] -t/2, xList[i+1][0] + t/2 ])

	if xList[-1][1] < lencontig - t :
		mList.append([xList[-1][1] -t/2, xList[-1][1] + t/2 ])

	return mList

def reportMisAssemblyIntervals2(xList, lencontig, t = 30, contigName=""):
	mList = []
	thres = 10
	gapAlignmentThres = 8000
	#gapAlignmentThres = 8000

	'''
	Production parameters initLenThres@ Helper, gapAlignment (8000, 12000) 
	i, precision, recall, TP_num, FP_num  : 3 	 0.084656 	 0.516129 	 16 	 173 	 0.145455 
	Quast analysis: (misassemblies, localMisassemblies) = (4, 6) 	,Original (18, 6)

	To find most TP, initThres, gapAlignment = (2000, 30)
	i, precision, recall, TP_num, FP_num  : 3 	 0.013188 	 0.741935 	 23 	 1721 	 0.025915 
	Quast analysis: (misassemblies, localMisassemblies) = (2, 0) 	,Original (18, 6) ... one is due to circular, another one is insertion missing

	Production parameters initLenThres@ Helper, gapAlignment (8000, 8000) 

	'''

	xList.sort()

	if xList[0][0] > t: 
		mList.append([xList[0][0], xList[0][0]+t/2]) 

	for i in range(len(xList) -1):
		d1, d2 = xList[i][1] - xList[i][0] , xList[i+1][1] - xList[i+1][0]
		if d1> thres and d2 > thres:
			diff = xList[i+1][0] - xList[i][1]
			if  diff < 0:
				mList.append([xList[i+1][0],xList[i][1] ])
			elif diff > gapAlignmentThres:
				mList.append([xList[i][1] - t/2, xList[i][1]])
				#mList.append([xList[i+1][0], xList[i+1][0] + t/2])


			#if 1000 < diff < 10000:
		#		print contigName, xList[i][1], xList[i+1][0], diff
			#else : 
		#		mList.append([xList[i][1] -t/2, xList[i+1][0] + t/2 ])
	if xList[-1][1] < lencontig - t :
		mList.append([xList[-1][1] -t/2, xList[-1][1]])

	return mList

def unitTest():
	xList = [[1,3], [4,7], [2,6]]
	intervalUnion(xList) 

	xList = [[1,3], [4,7]]
	intervalUnion(xList) 

	xList = [[1,3], [2,6]]
	intervalUnion(xList) 

	xList = [[1,3], [9,10], [2,6]]
	intervalUnion(xList) 
	
def unitTest2():
	if False:
		for i in range(4):
			print "i", i ," ============="
			xList = [[1,3], [4,7], [2,6], [1,2]]
			print intervalCover(xList, i) 

			xList = [[1,3], [4,7]]
			print intervalCover(xList, i) 

			xList = [[1,3], [2,6]]
			print intervalCover(xList, i) 

			xList = [[1,3], [9,10], [2,6]]
			print intervalCover(xList, i) 
			
			xList = [[1,3], [9,10], [2,6], [1.5, 4]]
			print intervalCover(xList, i) 

			xList = [[1,3], [9,10], [2,6], [1.5, 4]]
			coverList=  intervalCover(xList, i)
			print "coverList", coverList
			print "report", reportMisAssemblyIntervals(coverList, 10, i)
	else:
		xList = [[1, 10371], [10371, 18994], [19335, 33907], [27959, 33907]]
		coverList=  intervalCover(xList, 30)
		print "xList", xList
		print "coverList", coverList
		print "report", reportMisAssemblyIntervals(coverList, 33907, 30)

#unitTest2()

















