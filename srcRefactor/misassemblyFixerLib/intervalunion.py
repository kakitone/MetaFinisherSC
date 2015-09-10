
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

	#print xList, B
	return B

def unitTest():
	xList = [[1,3], [4,7], [2,6]]
	intervalUnion(xList) 

	xList = [[1,3], [4,7]]
	intervalUnion(xList) 

	xList = [[1,3], [2,6]]
	intervalUnion(xList) 

	xList = [[1,3], [9,10], [2,6]]
	intervalUnion(xList) 
	
#unitTest()
