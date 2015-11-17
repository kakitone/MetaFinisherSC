from itertools import groupby
from operator import itemgetter

import abunHouseKeeper

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import houseKeeper

import json

from multiprocessing import Pool
class seqGraphNodeWt(graphLib.seqGraphNode):
    def __init__(self, nodeIndex):
        graphLib.seqGraphNode.__init__(self, nodeIndex)
        self.nodeWt = 0 
            
        
    def getNextNodesFromIndices(self):
        returnList = []
        for eachitem in self.listOfNextNodes:
            returnList.append(eachitem[0])
        return returnList
    
    def getPrevNodesFromIndices(self):
        returnList = []
        for eachitem in self.listOfPrevNodes:
            returnList.append(eachitem[0])
        return returnList

        
class seqGraphWt(graphLib.seqGraph):
    def __init__(self, numberOfNodes):
        self.graphNodesList = [seqGraphNodeWt(i) for i in range(numberOfNodes)] 
    def formReportName(self, index):
        return str(index) + "_" +str(self.graphNodesList[index].nodeWt)

    def loadFromFile(self, folderName, fileName):
        
        f = open(folderName + fileName, 'r')
        numberOfNodes = 0 
        tmp = f.readline().rstrip()
        while (len(tmp) > 0):
            tmp = f.readline().rstrip()
            numberOfNodes += 1
        f.close()
        
        self.graphNodesList = [seqGraphNodeWt(i) for i in range(numberOfNodes)]
        
        f = open(folderName + fileName, 'r')
        tmp = f.readline().rstrip()
        runningIndex = 0
        while (len(tmp) > 0):
            dataList = tmp.split(';')
            

            for i in range(6):
                if len(dataList[i]) > 0:
                    myList = dataList[i].split(',')
                else:
                    myList = []

                    
                if i == 0:
                    self.graphNodesList[runningIndex].nodeIndex = int(dataList[i])
                elif i == 1:
                    self.graphNodesList[runningIndex].nodeIndexList = []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].nodeIndexList.append(int(eachitem))
                        
                elif i == 2:
                    self.graphNodesList[runningIndex].overlapList = []
                    for eachitem in myList:
                        self.graphNodesList[runningIndex].overlapList.append(int(eachitem))
                elif i == 3 or i == 4:
                    for eachitem in myList:
                        mydata = eachitem.split('-')
                        if i == 3 :
                            self.graphNodesList[runningIndex].listOfPrevNodes.append([int(mydata[0]), int(mydata[1])])
                        elif i == 4:
                            self.graphNodesList[runningIndex].listOfNextNodes.append([int(mydata[0]), int(mydata[1])])
                elif i == 5:
                    self.graphNodesList[runningIndex].visited = int(dataList[i])
            
                
            tmp = f.readline().rstrip()
            runningIndex = runningIndex + 1 
        f.close()
        
    def findConnectedComponents(self):
        for eachitem in self.graphNodesList:
            eachitem.visited = False 
        connectedList = []
        
        for eachnode in self.graphNodesList:
            if len(eachnode.nodeIndexList) > 0:
                tmpList = []
                
                if eachnode.visited == False:
                    stack = [eachnode]
                    tmpList = [eachnode.nodeIndex]
                    while len(stack) > 0:
                        currenttmp = stack.pop(0)
                        currenttmp.visited = True
                        
                        for eachsubIndex in currenttmp.listOfNextNodes:
                            eachsub = self.graphNodesList[eachsubIndex[0]]
                            if eachsub.visited == False:
                                stack.append(eachsub)
                                tmpList.append(eachsub.nodeIndex)
                    connectedList.append(tmpList)
                                
        
        return connectedList
    
        

class seqGraphDynamic(graphLib.seqGraph):
    
    def initAdv(self):
        self.N1 = len(self.graphNodesList)
        self.runningCtr = 0
        self.findAdjList()
        self.mapDummyToRealDic = {}
        self.xResolvedSimplifiedList = []
        
    def findAdjList(self):
        self.adj = [[] for i in range(self.N1)]
        
        for i in range(self.N1):
            for eachnext in self.graphNodesList[i].listOfNextNodes: 
                if len(self.graphNodesList[eachnext[0]].nodeIndexList) > 0:
                    self.adj[i].append(eachnext[0])
            
    def clearOut(self, u):
        outList  = []

        for eachnext in self.graphNodesList[u].listOfNextNodes:
            outList.append(eachnext[0])
        
        for x in outList:
            self.removeEdge(u,x)
        
    def clearIn(self, v ):
        inList  = []

        for eachprev in self.graphNodesList[v].listOfPrevNodes:
            inList.append(eachprev[0])
        
        for x in inList:
            self.removeEdge(x, v)

    def clearVisitStatus(self):
        for eachnode in self.graphNodesList:
            eachnode.visited = False
        
    def reachable(self, i, j):
        self.clearVisitStatus()
        depthMax = 3 
        
        queue = [i]
        levelQ = [0] 
        
        self.graphNodesList[i].visited = True
        
        canReach = False
        
        while len(queue) > 0:
            current = queue.pop(0)
            currentLevel = levelQ.pop(0)
            
            if current == j :
                canReach = True
                break
            
            if currentLevel >= depthMax:
                break
            
            
            for eachnext in self.graphNodesList[current].listOfNextNodes:
                index = eachnext[0]
                if self.graphNodesList[index].visited == False:
                    queue.append(index)
                    self.graphNodesList[index].visited = True
                    levelQ.append(currentLevel+1)
                    
        self.clearVisitStatus()
        return canReach
    
    def logEdges(self, folderName, stagename):
        print "Logging edges"
        logList = []
        
        if stagename == "XResolution":
            mapDummyToRealDic =self.readInJSON(folderName, "mapDummyToRealDic.json")
        else:
            mapDummyToRealDic = {}

        for eachnode in self.graphNodesList:
            tmpNodeIndexList = []
            for kk in eachnode.nodeIndexList:
                if kk >= self.N1:
                    tmpNodeIndexList += mapDummyToRealDic[str(kk-self.N1)][1]
                else:
                    tmpNodeIndexList += [kk]

            if len(tmpNodeIndexList) >= 2:
                for i in range(len(tmpNodeIndexList)-1):
                    currentName = tmpNodeIndexList[i]
                    nextName =  tmpNodeIndexList[i+1]
                    cName =  abunHouseKeeper.parseIDToName(currentName,'C',0)
                    nName =  abunHouseKeeper.parseIDToName(nextName,'C',0)
                    logList.append([cName, nName])

        with open( folderName + stagename + ".json", 'w') as f:
            json.dump(logList, f)    

    def condenseEdgeRemove(self, G_ContigRead, folderName, mummerLink, contigFilename):
        print "condenseEdgeRemove"
        thresPass = 100
        thresForStrangeCut = 5000
        ### kkdebug

        toRemoveList = []
        
        for eachnode in self.graphNodesList:
            if len(eachnode.nodeIndexList) > 0:
                if len(eachnode.listOfNextNodes) ==1  :
                    nextNodeIndex = eachnode.listOfNextNodes[0][0]
                    nextNode= self.graphNodesList[nextNodeIndex]
                    if len(nextNode.listOfPrevNodes) == 1 : 
                        currentName = eachnode.nodeIndex
                        nextName =  nextNode.nodeIndex

                        contigReadPaths = findAllPathK(currentName,nextName, G_ContigRead, 5)

                        cName =  abunHouseKeeper.parseIDToName(currentName,'C',0)
                        nName =  abunHouseKeeper.parseIDToName(nextName,'C',0)

                        noGoNext = self.readInJSON(folderName, "noGoNext.json")
                        noGoPrev = self.readInJSON(folderName, "noGoPrev.json")

                        overlap = [-1, -1]
                        ctr = 0 

                        for eachpath in contigReadPaths:
                            if len(eachpath) > 2: 
                                ctr = ctr + 1 
                                
                            elif len(eachpath) == 2:     
                                
                                contigName = cName
                                leftSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)

                                contigName = nName
                                rightSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)
                                
                                overlap = IORobot.align(leftSeg, rightSeg, folderName, mummerLink)


                        if ctr <= thresPass and  (cName in noGoNext or nName in noGoPrev or overlap[0] > thresForStrangeCut ):
                    
                            self.removeEdge(currentName, nextName)
                            toRemoveList.append([currentName, nextName])


        ### kkdebug
        #with open( "dataFolder/toRemoveList.json", 'w') as f:
        #    json.dump(toRemoveList, f)    

        self.findAdjList()

    def readInJSON(self, folderName, filename):

        json_data = open(folderName + filename, 'r')
        dataItem = json.load(json_data)
        return dataItem
    
    def transitiveReduction(self,folderName, mummerPath, contigFile, readFile, G_ContigRead):
        for i in range(self.N1):
            for j in self.adj[i]:
                self.removeEdge(i, j)
                
                canReach = self.reachable(i,j)
                
                if canReach == False:
                    self.insertEdge(i, j, 1997)

                if False:
                    contigPaths = findAllPathK(i,j, self,3)
                    contigReadPaths = findAllPathK(i, j ,G_ContigRead ,5)

                    directPathList = []
                    for eachp in contigReadPaths:
                        ck = False
                        for pitem in eachp[1:-1] :
                            if pitem < self.N1:
                                ck = True
                        if ck == False:
                            directPathList.append(eachp)

                    indirectPathList = []
                    for eachp in contigPaths:
                        if len(eachp) > 2:
                            tmpPath = [eachp[0]]
                            for pp in range(len(eachp) - 1):
                                startI, endI = eachp[pp], eachp[pp+1]
                                pathList = findAllPathK(startI, endI , G_ContigRead,3)
                                tmpPath = tmpPath + pathList[0][1:]
                            
                            indirectPathList.append(tmpPath)

                    formPathSeq(folderName, mummerPath, directPathList, indirectPathList, contigFile, readFile)

                    toDelete = decideCut(folderName, mummerPath)
                    if not toDelete : 
                        self.insertEdge(i,j, 1997)
        self.findAdjList()



    def doubleEdgeReduction(self):
        for u in range(self.N1):
            for v in self.adj[u]:
                if u in self.adj[v]:
                    self.removeEdge(u,v)
                    self.removeEdge(v,u)


        self.findAdjList()


    def bipartiteResolve(self, resolvedList):
        for e in resolvedList:
            u, v =e[0], e[-1]
            self.clearOut(u)
            self.clearIn(v)
            
            self.insertEdge(u,v,1997)
    

    def bipartiteLocalResolve(self, resolvedList, inList, outList, folderName):

        #noGoNext = self.readInJSON(folderName, "noGoNext.json")
        #noGoPrev = self.readInJSON(folderName, "noGoPrev.json")

        if len(resolvedList) > 0:
            for u in inList:
                self.clearOut(u/2)

            for v in outList:
                self.clearIn(v/2)

            for e in resolvedList:
                u, v =e[0], e[-1]
                cName =  abunHouseKeeper.parseIDToName(u,'C',0)
                nName =  abunHouseKeeper.parseIDToName(v,'C',0)

                #if not cName in noGoNext and not nName in noGoPrev:  
                self.insertEdge(u,v,1997)


    def symGraph(self):
        for u in range(len(self.graphNodesList)):
            uComp = findMyComp(u)
            
            for vI in self.graphNodesList[u].listOfNextNodes:
                v = vI[0]
                vComp = findMyComp(v)

                if not graphLib.nameInEdgeList(vComp, self.graphNodesList[uComp].listOfPrevNodes):
                    self.removeEdge(u,v)

            for vI in self.graphNodesList[u].listOfPrevNodes:
                v = vI[0]
                vComp = findMyComp(v)

                if not graphLib.nameInEdgeList(vComp, self.graphNodesList[uComp].listOfNextNodes):
                    self.removeEdge(v,u)




    def xResolve(self, xResolvedList):
        self.mapDummyToRealDic = {}
        self.xResolvedSimplifiedList = []
        
        for i in range(self.N1):
            for y in xResolvedList[i]:
                
                inNode = y[0]
                outNode = y[-1]
                
                for eachprev in self.graphNodesList[i].listOfPrevNodes:
                    key =  eachprev[0]
                    if key >= self.N1 and self.mapDummyToRealDic[key-self.N1][0] == y[0]:
                        inNode = key 
                
                for eachnext in self.graphNodesList[i].listOfNextNodes:
                    key =  eachnext[0]
                    if key >= self.N1 and self.mapDummyToRealDic[key-self.N1][0] == y[-1]:
                        outNode = key 
                
                
                if len(self.graphNodesList[i].listOfPrevNodes) == 1 and len(self.graphNodesList[i].listOfNextNodes) ==1: 

                    self.xResolvedSimplifiedList.append([inNode, i])
                    self.xResolvedSimplifiedList.append([i, outNode])
                else:
                    self.graphNodesList.append(graphLib.seqGraphNode(self.N1+self.runningCtr))
                    
                    if False:
                        toRemoveList = []

                        for outgo in self.graphNodesList[inNode].listOfNextNodes:
                            toRemoveList.append([inNode, outgo[0]])    

                        for income in self.graphNodesList[outNode].listOfPrevNodes:    
                            toRemoveList.append([income[0], outNode])

                        for eachpair in toRemoveList:
                            self.removeEdge(eachpair[0], eachpair[1]) 
                    else:
                        self.removeEdge(inNode,i )
                        self.removeEdge(i, outNode)

                    self.insertEdge(inNode, self.N1 + self.runningCtr, 1997)
                    self.insertEdge(self.N1 + self.runningCtr, outNode, 1997)

                    self.mapDummyToRealDic[self.runningCtr] = [ [i] , self.graphNodesList[i].nodeIndexList ]
                    
                    self.xResolvedSimplifiedList.append([inNode, self.N1 + self.runningCtr])
                    self.xResolvedSimplifiedList.append([self.N1 + self.runningCtr, outNode])
                    self.runningCtr = self.runningCtr + 1
                

def findMyComp(u):

    if u %2 == 0:
        myComp = u+1
    else:
        myComp = u-1
    return myComp
    
def findAllReachable(i, N1, G):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        
    myList = []
    tmpList = [G.graphNodesList[i]]
    neighborList = []
    while len(tmpList) > 0:
        current = tmpList.pop(0)
        for eachChildIndex in current.listOfNextNodes:
            eachChild = G.graphNodesList[eachChildIndex[0]]
            
            if eachChild.nodeIndex >= N1:
                if eachChild.visited == False:
                    eachChild.visited = True
                    tmpList.append(eachChild)
            else:
                neighborList.append(eachChild)
    
    connectedNodeIndexList = []
    
    for eachitem in neighborList:
        connectedNodeIndexList.append(eachitem.nodeIndex)
    
    connectedNodeIndexList.sort()

    for key, items in groupby(connectedNodeIndexList):
        myList.append(key)
    
    return myList 



def filterEdge(adjacencyList, folderName, contigFilename):
    lenDic = IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")
    thresFoPhase = 2000
    smallList, largeList = [], []
    for eachitem in lenDic:
        id = abunHouseKeeper.parseEdgeNameToID(eachitem, 'C')
        if lenDic[eachitem] < thresFoPhase:
            smallList.append(id)
        else:
            largeList.append(id)
    
    newAdjacencyList = [[] for i in range(len(adjacencyList))]
    
    for i in largeList:
        for eachitem in adjacencyList[i]:
######## IMPORTANT:
            if  eachitem in largeList and eachitem / 2 != i / 2:
######## NEED TO REMOVE IN PRODUCTION if True
                newAdjacencyList[i].append(eachitem)
    
    
    print "len(smallList)  , len(largeList): ", len(smallList)  , len(largeList)
    print "lenDic: ", lenDic
    
    for eachitem in newAdjacencyList:
        print "newAdjacencyList :", eachitem 
        
    return newAdjacencyList
        

def addIndicesToReachable(x, G, N1):
    
    
    tmp = int(x.split('_')[0])

    findAllReachable(tmp , N1, G)
    
    for eachnode in G.graphNodesList:
        if eachnode.visited == True and eachnode.nodeIndex >= N1:
            eachnode.visitLabelList.append(x)


def markReachableIndices(G, Grev, rIn, rOut, N1):
    print "markReachableIndices"
    for eachitem in G.graphNodesList:
        eachitem.visitLabelList = []
        eachitem.visited = False
    
    for eachitem in Grev.graphNodesList:
        eachitem.visitLabelList = []
        eachitem.visited = False
    
    for x in rIn:
        addIndicesToReachable(x, G, N1)
    
    # # Bug : wrong indexing 

    for y in rOut:
        addIndicesToReachable(y, Grev, N1)
        
    for i in range(len(Grev.graphNodesList)):
        if len(Grev.graphNodesList[i].visitLabelList) > 0:
            G.graphNodesList[i].visitLabelList = G.graphNodesList[i].visitLabelList + Grev.graphNodesList[i].visitLabelList  

def formReverseGraph(G):
    nNode = len(G.graphNodesList)
    Grev = seqGraphWt(nNode)
    for i in range(nNode):
        for j in range(nNode):
            haveInserted = graphLib.nameInEdgeList(j, G.graphNodesList[i].listOfNextNodes) 
            if haveInserted:      
                Grev.insertEdge(j, i, 100)
    return Grev

def formReverseGraphFast(G):
    nNode = len(G.graphNodesList)
    Grev = seqGraphWt(nNode)
    for i in range(nNode):
        tmpList = G.graphNodesList[i].listOfNextNodes
        Grev.graphNodesList[i].listOfPrevNodes = tmpList
        
        tmpList = G.graphNodesList[i].listOfPrevNodes
        Grev.graphNodesList[i].listOfNextNodes = tmpList


    return Grev


def markInsideNodes(G, kkIn, kkOut):

    print "markInsideNodes"
    singleMissList = [[int(kkIn[0].split('_')[0])] , [int(kkIn[1].split('_')[0])], [int(kkOut[0].split('_')[0])], [int(kkOut[1].split('_')[0])]]
    allPassList = []
    
    for eachitem in G.graphNodesList:
        result = checkMiss(eachitem, kkIn, kkOut)
        if result == -1: 
            allPassList.append(eachitem.nodeIndex)
        elif result == 0:
            singleMissList[1].append(eachitem.nodeIndex)
        elif result == 1:
            singleMissList[0].append(eachitem.nodeIndex)
        elif result == 2:
            singleMissList[3].append(eachitem.nodeIndex)
        elif result == 3:
            singleMissList[2].append(eachitem.nodeIndex)
        
        # elif 0 <= result <= 3:
        #    singleMissList[result].append(eachitem.nodeIndex)
        
    return singleMissList, allPassList


def checkMiss(myNode, rIn, rOut):
    index = 1997
    myList = myNode.visitLabelList

    combineCheckList = rIn + rOut
    if set(myList) == set(combineCheckList):
        return -1
    for i in range(4):
        if set(myList) == set(combineCheckList[0:i] + combineCheckList[i + 1:]):
            return i
            break
    return index

 
def markStartEndNodes(G, rIn, rOut, singleMissList, allPassList):
    print "markStartEndNodes"
    myStartIndex, myEndIndex = -1, -1
    
    print singleMissList 
    print allPassList
    for i in allPassList:
        nextList = G.graphNodesList[i].getNextNodesFromIndices()
        prevList = G.graphNodesList[i].getPrevNodesFromIndices()
        
        counter = 0 
        for j in range(2):
            if len(set(prevList).intersection(set(singleMissList[j]))) > 0:
                counter = counter + 1
        
        
        for j in range(2, 4):
            if len(set(nextList).intersection(set(singleMissList[j]))) > 0:
                counter = counter + 10
        
        if counter / 10 == 2 and counter % 10 == 0:
            myEndIndex = i
        
        if counter / 10 == 0 and counter % 10 == 2:
            myStartIndex = i
            
        if counter == 22:
            myStartIndex = i
            myEndIndex = i
            break
    
    return myStartIndex, myEndIndex
            

def BFS(x, y, G, N1):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        eachnode.backPtr = -1
        
    tmpList = [G.graphNodesList[x]]
    starter = True
    while len(tmpList) > 0:
        current = tmpList.pop(0)
        if current.nodeIndex >= N1 or starter:
            starter = False
            for i in current.getNextNodesFromIndices():
                if G.graphNodesList[i].visited == False:
                    G.graphNodesList[i].visited = True
                    G.graphNodesList[i].backPtr = current.nodeIndex
                    tmpList.append(G.graphNodesList[i])
            
        
    # print "end BFS", y ,x
    z = y
    tmpList.append(y)           
    while z != x:
        z = G.graphNodesList[z].backPtr
        tmpList.append(z)
    
    # print "end Back Ptr"

    return tmpList[::-1]

def markInterior(G , myStartIndex, myEndIndex, N1):
    print "markInterior"
    myPathway = BFS(myStartIndex, myEndIndex , G, N1)
    return myPathway

   
def markFlankingRegion(G, rIn, rOut, myStartIndex, myEndIndex, N1):
    
    print "markFlankingRegion"
    print myStartIndex, myEndIndex
    myPathwayList = [[] for i in range(4)]
    
    for i in range(2):
        myPathwayList[i] = BFS(rIn[i], myStartIndex, G, N1)
        myPathwayList[i + 2] = BFS(myEndIndex, rOut[i], G, N1)
        
    return myPathwayList

def DFS(G, x, N1, startIndex, endIndex, mypath):
    
    if (x >= N1 or x == startIndex ) and len(mypath) < 3:
        if G.graphNodesList[x].visited == False:
            G.graphNodesList[x].visited = True
            for eachChild in G.graphNodesList[x].listOfNextNodes:
                #print "startIndex, endIndex, eachChild[0], x", startIndex, endIndex, eachChild[0], x
                if G.graphNodesList[eachChild[0]].visited == False:
                    returnpath = DFS(G, eachChild[0], N1, startIndex, endIndex, mypath + [x])
                else:
                    returnpath = None

                if returnpath != None :
                    return returnpath
                
        return None
            
    elif x == endIndex:
        print mypath + [x]
        return mypath + [x]
    else:
        return None
        
        
def recCheck(G, i , N1, counter, contigList):
    if counter > 0:
        for eachitem in G.graphNodesList[i].listOfNextNodes:
            if G.graphNodesList[eachitem[0]].nodeIndex < N1:
                contigList.append(eachitem[0])
                continue
            recCheck(G, eachitem[0] , N1, counter - 1, contigList)
            
def checkFourHoppers(G, i , N1):
    contigList = []
    recCheck(G, i , N1, 3, contigList)
    
    returnList = []
    contigList.sort()
    for key, items in groupby(contigList):
        returnList.append(key)
    return returnList
                        
def debugGraphPath(startIndex, endIndex, G, N1):
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        
    print startIndex , checkFourHoppers(G, startIndex , N1)

    
def markAssociatedReads(G, singleMissList, allPassList):
    print "markAssociatedReads"
    flankingList = [[] for i in range(4)]
    repeatList = []
    for eachitem in G.graphNodesList:
        nextList = eachitem.getNextNodesFromIndices()
        prevList = eachitem.getPrevNodesFromIndices()
        
        
        if len(set(nextList).intersection(set(allPassList))) > 0 :
            nextInt = True
        else:
            nextInt = False
        
        if len(set(prevList).intersection(set(allPassList))) > 0 :
            prevInt = True
        else:
            prevInt = False
        
        
        missArr = [False for i in range(4)]
        for j in range(2):
            if len(set(prevList).intersection(set(singleMissList[j]))) > 0 :
                missArr[j] = True
            if len(set(nextList).intersection(set(singleMissList[j + 2]))) > 0:
                missArr[j + 2] = True
        
        for j in range(2):
            if missArr[j] and prevInt:
                flankingList[j].append(eachitem.nodeIndex)
            elif missArr[j + 2] and nextInt:
                flankingList[j + 2].append(eachitem.nodeIndex)
            elif prevInt and nextInt:
                repeatList.append(eachitem.nodeIndex)
        
    return flankingList, repeatList
            

def checkPathLength(path, G, N1, folderName):
    
    lenDicRR = IORobot.obtainLength(folderName, "phasingSeedName_Double.fasta")
    sumLength = 0
    overlapLength = 0
    for index, i in zip(path, range(len(path))):
        header = "Read" + str((index - N1) / 2) + "_"
        if (index - N1) % 2 == 0:
            header = header + "p"
        else:
            header = header + "d"
        print "lenDicRR[header], ", lenDicRR[header], header 
        print (index - N1) * 2 + 1, (index - N1) * 2 + 2
        sumLength = sumLength + lenDicRR[header]
        
        if i != len(path) - 1:
            for eachnext in G.graphNodesList[index].listOfNextNodes:
                if eachnext[0] == path[i + 1]:
                    overlapLength = overlapLength + eachnext[1]
                    break 
    print sumLength, overlapLength, sumLength - overlapLength

def findPathBtwEnds(folderName, leftCtgIndex, rightCtgIndex, contigReadGraph, N1):
    
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    
    return findPathBtwEndsFast(folderName, leftCtgIndex, rightCtgIndex, G, N1)

def findPathBtwEndsFast(folderName, leftCtgIndex, rightCtgIndex, G, N1):
    readList = []
    
    for eachnode in G.graphNodesList:
        eachnode.visited = False
        
    x, startIndex, endIndex, mypath = leftCtgIndex,leftCtgIndex ,rightCtgIndex, []
    returnPath = DFS(G, x, N1, startIndex, endIndex, mypath)
    
    if returnPath == None:
        return None
    else:

        if returnPath[0] == startIndex:
            returnPath.pop(0)
        
        if returnPath[-1] == endIndex:
            returnPath.pop(-1)
        
        readList = returnPath 
        
        return readList 


def findAllPathK(u,v,G,k):
    allPathList = BFS_revisit(u,v,G,k)
    return allPathList

def BFS_revisit(u, v, G, k):
    resultList = []
    isEdge = graphLib.nameInEdgeList(v, G.graphNodesList[u].listOfNextNodes)

    if isEdge:
        resultList.append([u,v])

    if k>1:
        for w in G.graphNodesList[u].listOfNextNodes:
            returnList = BFS_revisit(w[0],v,G,k-1)
            for p in returnList:
                tmpList = [u]
                for i in p:
                    tmpList.append(i)

                resultList.append(tmpList)

    return resultList

def formPathSeq(folderName, mummerPath, directPathList, indirectPathList, contigFile, readFile):
    '''
    Input : directPathList, indirectPathList, contigFile, readFile
    Output: directPath.fasta, indirectPath.fasta
    '''

    contigList = IORobot.readContigsFromFile(folderName,contigFile)
    readList = IORobot.readContigsFromFile(folderName,readFile)
    
    directPathSeqList =  IORobot.pathListToSeqListTransform(directPathList, contigList, readList, mummerPath, folderName)    
    indirectPathSeqList =  IORobot.pathListToSeqListTransform(indirectPathList, contigList, readList, mummerPath, folderName)    

    IORobot.writeSegOut(directPathSeqList,folderName,"directPath.fasta")
    IORobot.writeSegOut(indirectPathSeqList,folderName,"indirectPath.fasta")

def decideCut(folderName, mummerPath):
    
    '''
    Input : directPath.fasta, indirectPath.fasta
    Output : toDelete 
    '''
    thres = 50
    
    if True:
        alignerRobot.useMummerAlign(mummerPath, folderName, \
            "indirectvsdirect", "indirectPath.fasta", "directPath.fasta", specialForRaw = False, specialName = "", refinedVersion= True)
    
    dataList =  alignerRobot.extractMumData(folderName , "indirectvsdirectOut")
    lenDic = IORobot.obtainLength(folderName, "directPath.fasta")

    ctr =0 
    ctrindirect = 0 

    dataList.sort(key = itemgetter(-1))

    toDelete = True

    for key, items in groupby(dataList, itemgetter(-1)):
        print "key", key 
        ctr = ctr + 1
        isFound = False
        for eachitem in items:
            if eachitem[2] < thres and eachitem[3] > lenDic[key] - thres:
                isFound = True

        if isFound:
            ctrindirect = ctrindirect + 1


    epsilon = 1.1

    print "ctrindirect, ctr", ctrindirect, ctr

    if ctrindirect*1.0/ctr < (1- epsilon):
        toDelete = False
    else:
        toDelete = True


    return toDelete

def parallelGapLookUp(resolvedList,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename):
    if houseKeeper.globalRunMPI == False:
        p = Pool(processes=houseKeeper.globalParallel)
        results = []
        
        for eachmatchpair in resolvedList:
            results.append(p.apply_async(singleGapLookUp, args=(eachmatchpair,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename)))

        outputlist = [itemkk.get() for itemkk in results]
        print  len(outputlist)
        p.close()
    else:
        from mpi4py import MPI
        from mpi4py.MPI import ANY_SOURCE

        comm = MPI.COMM_WORLD
        me = comm.Get_rank()
        numberOfWorkers = comm.Get_size() - 1

        results = []
        outputlist = []

        for eachmatchpair in resolvedList:
            results.append([eachmatchpair,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename])

        for i in range(len(results)):
            data = results[i]
            data.insert(0, "gapjob")
            print "master sender", data[0] 
            comm.send(data, dest=(i%numberOfWorkers) +1)
            

        for i in range(len(results)):    
            data = comm.recv(source=ANY_SOURCE)
            print "master receiver", data[0:3]
            outputlist.append(data)



    return outputlist

def singleGapLookUp(eachmatchpair,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename):
    #print eachmatchpair
    leftCtgIndex ,rightCtgIndex, leftEnd, rightStart, middleContent = eachmatchpair[0],eachmatchpair[-1],0,0,""
    
    succReadsList = []
    G = seqGraphWt(0)
    G.loadFromFile(folderName, contigReadGraph)
    succReadsList = BFS(leftCtgIndex,rightCtgIndex, G, N1)

    if len(succReadsList) > 0:
        succReadsList.pop(0)
        succReadsList.pop(-1)
    else:
        print "interesting item for future study"

    print "succReadsList" , succReadsList
    
    if len(succReadsList) == 0:
        contigName = abunHouseKeeper.parseIDToName(leftCtgIndex, 'C', N1)
        leftSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)

        contigName = abunHouseKeeper.parseIDToName(rightCtgIndex, 'C', N1)
        rightSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)
        
        overlap = IORobot.alignWithName(leftSeg, rightSeg, folderName, mummerLink, str(leftCtgIndex) + "_" + str(rightCtgIndex) )
        
        print "overlap contig : ", overlap
        
        leftEnd = len(leftSeg) - overlap[0]
        middleContent = ""
        
    else:
        
        contigName = abunHouseKeeper.parseIDToName(leftCtgIndex, 'C', N1)
        print contigName
        leftSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)
        
        readName = abunHouseKeeper.parseIDToName(succReadsList[0], 'R', N1)
        print readName
        rightSeg  = IORobot.myRead(folderName, readsetFilename + "_Double.fasta", readName)
        
        overlap = IORobot.alignWithName(leftSeg, rightSeg, folderName, mummerLink, str(leftCtgIndex) + "_" + str(rightCtgIndex) )
        
        print "overlap start read : ", overlap
        
        leftEnd = len(leftSeg) - overlap[0]
        
        middleContent = ""
        
        for i in range(len(succReadsList)-1):
            readName = abunHouseKeeper.parseIDToName(succReadsList[i], 'R', N1)
            leftSeg  = IORobot.myRead(folderName, readsetFilename + "_Double.fasta", readName)
        
            readName = abunHouseKeeper.parseIDToName(succReadsList[i+1], 'R', N1)
            rightSeg  = IORobot.myRead(folderName, readsetFilename + "_Double.fasta", readName)
            
            overlap = IORobot.alignWithName(leftSeg, rightSeg, folderName, mummerLink, str(leftCtgIndex) + "_" + str(rightCtgIndex) )
            print "overlap middle read : ", overlap
            middleContent = middleContent + leftSeg[0:len(leftSeg)-overlap[0]] 
        
        
        readName = abunHouseKeeper.parseIDToName(succReadsList[-1], 'R', N1)
        leftSeg  = IORobot.myRead(folderName, readsetFilename + "_Double.fasta", readName)
        
        contigName = abunHouseKeeper.parseIDToName(rightCtgIndex, 'C', N1)
        rightSeg = IORobot.myRead(folderName, contigFilename + "_Double.fasta", contigName)
        
        overlap = IORobot.alignWithName(leftSeg, rightSeg, folderName, mummerLink, str(leftCtgIndex) + "_" + str(rightCtgIndex) )
        print "overlap end read : ", overlap
        
        middleContent = middleContent + leftSeg[0:len(leftSeg)-overlap[0]]

    return [leftCtgIndex ,rightCtgIndex, leftEnd, rightStart, middleContent]
    

