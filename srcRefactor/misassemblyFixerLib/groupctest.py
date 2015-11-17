import numpy as np
import random
import scipy.stats
'''
Simulation code : 

1. Data generation : 
	Input : n=100, 
	Output : (x1,x2,x3,...,x100), (g1,g2,...,g100)

2. Define decision boundary : 
	Input = (g1,g2,g3,...,g100)
	Output = (b1,b2,...,b99)

3. Algorithms 1,2,3: 
	Input: (x1,x2,x3,...,x100)
	Output: (p1,p2,...,p99)

4. Calcuate precision and recall 
	Input : b, p
	Output : precision, recall
'''

def dataGeneration():
	n = 100
	m = [20, 50]
	x = [0 for i in range(n)]
	g = [0 for i in range(n)]

	x[0] = np.random.poisson(m[0])
	g[0] = 0
	for i in range(1, n):
		rand = random.random()
		if rand > 0.1:
			prev = g[i-1] 
			x[i] = np.random.poisson(m[prev])
			g[i] = prev 
		else:
			notprev = 1 if g[i-1] == 0 else 0
			x[i] = np.random.poisson(m[notprev])
			g[i] = notprev

	return x, g

def findDecisionBoundaryFromGT(g):
	b = [-1 for i in range(len(g)-1)]

	for i in range(len(g)-1):
		if g[i] == g[i+1]:
			b[i] = 0
		else:
			b[i] = 1
	return b

def predictDecisionBoundary(x, method="basic"):
	p = []
	if method=="basic":
		p = [1 for i in range(len(x)-1)]
	elif method == "succCTest":
		p = succCTest(x)
	elif method == "groupCTest":
		p = groupCTest(x)

	return p

def succCTest(x):
	# scipy.stats.binom_test(x, n=None, p=0.5)
	alpha = 0.05
	p = [1 for i in range(len(x)-1)]
	for i in range(len(x) -1):
		pvalue = scipy.stats.binom_test(x[i][0], n=x[i][0]+x[i+1][0], p=0.5)
		p[i] = 1 if pvalue < alpha else 0

	return p

def groupCTest(x):
	alpha = 0.05
	p = [1 for i in range(len(x) -1) ]
	u = [0]

	count = 0
	while (np.sum(u) != len(u)):
		#print count
		count += 1
		gp = findUpdateGrouping(p, x)
		#print len(gp), len(x)
		u = performGpCTest(gp, alpha)
		#print len(u)
		p = updateCurrentVector(p, u)

	return p 

def findUpdateGrouping(p, x):
	gp = [] 
	# x = [10, 20, 21, 19, 20, 26]
	# p = [  1,  0 ,  0,  0,  1]
	# gp = [[10, 1], [80, 4], [26, 1]]

	prev, current = 0, 1



 	for i in range(len(p)):
 		if p[i] == 1:
			tmpsum, lensum = calculateSum(x[prev:current])
 			gp.append([tmpsum, lensum])
			prev = current
			current = prev + 1 
		elif p[i] == 0:
			current = current + 1

	if p[-1] == 0 :
		tmpsum, lensum = calculateSum(x[prev:])
		gp.append([tmpsum, lensum])
	else:
		gp.append([x[-1][0], x[-1][1]])

	return gp

def calculateSum(xlist):
	tmpsum = 0
	lensum = 0
	maxVal = 10**6 
	for eachitem in xlist:
		if tmpsum < maxVal:
			tmpsum = tmpsum + eachitem[0]
			lensum = lensum + eachitem[1]
		else:
			break
	return tmpsum, lensum

def performGpCTest(gp, alpha):
	u = [1 for i in range(len(gp)-1)] 
	# gp = [[10, 1], [80, 4], [26, 1]]
	for i in range(len(u)):
		p_tmp = gp[i][1]*1.0 / (gp[i][1] + gp[i+1][1])
		x1_tmp = gp[i][0] 
		n_tmp = gp[i][0] + gp[i+1][0]

		pvalue = scipy.stats.binom_test(x1_tmp, n=n_tmp, p=p_tmp)

		u[i] = 1 if pvalue < alpha else 0

	return u

def updateCurrentVector(p, u):
	pnew = [p[i] for i in range(len(p))]

	curr = 0
	for i in range(len(p)):
		if p[i] == 0:
			pnew[i] == 0
		elif p[i] == 1:
			pnew[i] = u[curr]
			curr = curr + 1

	return pnew

def calculatePrecisionAndRecall(b,p):
	T, P, TP = 0,0,0
	#assert(len(b) == len(p))
	for i in range(len(b)):
		if b[i] == 1:
			T += 1

		if p[i] == 1:
			P += 1

		if b[i] == 1 and p[i] == 1:
			TP += 1

	if  P == 0 :
		precision = 1 
	else:
		precision = TP*1.0/P

	if T ==0:
		recall = 1 
	else:
		recall =  TP*1.0/T
	
	return precision, recall

def mainFlow():
	numberOfRounds = 100
	methodList = ["basic", "succCTest", "groupCTest"]
	#methodList = ["groupCTest"]
	for method in methodList :
		precisionList, recallList = [], []
		for i in range(numberOfRounds):
			x, g = dataGeneration()
			b = findDecisionBoundaryFromGT(g)
			xnew = transformwithwtinfo(x)

			p = predictDecisionBoundary(xnew, method)
			precision, recall= calculatePrecisionAndRecall(b,p)
			precisionList.append(precision)
			recallList.append(recall)

		f1score = 2*np.mean(precisionList)*np.mean(recallList) /(np.mean(precisionList) +  np.mean(recallList))
		sdprec, sdrecall = np.std(precisionList) , np.std(recallList)
		print method + "-- precision, recall, std(prec), std(recall), F1score: %f , %f , %f, %f , %f"%(np.mean(precisionList), np.mean(recallList),sdprec, sdrecall, f1score)

def transformwithwtinfo(x):
	xnew = []
	for eachitem in x:
		xnew.append([eachitem, 1])
	return xnew 

#mainFlow()







