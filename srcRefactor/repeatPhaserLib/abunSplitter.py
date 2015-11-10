import json
import numpy as np
import os
from itertools import groupby
from operator import itemgetter
import copy
import sys
from random import shuffle
from multiprocessing import Pool

from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import IORobot
from finisherSCCoreLib import houseKeeper
from finisherSCCoreLib import nonRedundantResolver 

import associatedReadFinder
import readContigGraphFormer
import repeatFinder
import abunHouseKeeper
import abunGraphLib
import emalgo


### Abundance split and its subroutines
def obtainNonEmpty(repeatPairs):
    newRepeatPairs = []
    for eachitem in repeatPairs:
        if len(eachitem[0]) > 0 and len(eachitem[1]) > 0:
            newRepeatPairs.append(eachitem)
    
    return newRepeatPairs

def convert4to1base(i):
    return i/4

def convert4to2base(i):
    return i/2

def getCt(inList ,myCountDic):
    newInList = []
    for i in inList: 
        tmp1 = convert4to1base(i)
        tmp2 = convert4to2base(i)
        newInList.append([tmp2 , myCountDic["Segkk"+str(tmp1)]])
    return newInList 

def getCtAgg(inList ,myCountDic, Gnew, lenDic):
    newInList = []
    print "New getCtAgg", len(inList)
    for i in inList: 
        tmp1 = convert4to1base(i)
        tmp2 = convert4to2base(i)
        covTmp = 0
        lenTmp = 0

        print "len(Gnew.graphNodesList[tmp2].nodeIndexList), len(inList)", len(Gnew.graphNodesList[tmp2].nodeIndexList), len(inList)

        for eachindex in Gnew.graphNodesList[tmp2].nodeIndexList:
            
            if lenDic[abunHouseKeeper.parseIDToName(eachindex, 'C', len(lenDic))] > lenTmp :
                lenTmp = lenDic[abunHouseKeeper.parseIDToName(eachindex, 'C', len(lenDic))]
            
                name = "Segkk"+str(eachindex/2)
                covTmp = myCountDic[name]

        newInList.append([tmp2 , covTmp])

    return newInList 

def getCtTwoToOne(inList ,myCountDic):
    newInList = []
    for i in inList: 
        newInList.append([i , myCountDic["Segkk"+str(i/2)]])
    return newInList 

def getCtTwoToOneAgg(inList ,myCountDic, Gnew, lenDic):

    newInList = []
    for eachitem in inList:
        newInList.append(eachitem*2)

    return getCtAgg(newInList ,myCountDic, Gnew, lenDic)

def satisfyMatch(initem ,newOutList, sd):
    print "AbunLower, AbunUpper", abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLower,abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpper

    found = -1
    inIndex , inCt = initem[0], initem[1]
    # First check pt
    targetItem = -1
    for outitem in newOutList:
        outIndex, outCt = outitem[0] , outitem[1]
        print abs(outCt - inCt)
        if abs(outCt - inCt) < abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLower*sd and inIndex != outIndex:
            targetItem = outIndex
            
    # Second check pt
    
    rejection = False
    for outitem in newOutList:
        outIndex, outCt = outitem[0] , outitem[1]
        if outIndex != targetItem :
            if abs(outCt - inCt) <=abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpper*sd :
                rejection = True
                print "rejection ", abs(outCt - inCt)
                
    # Combined check 
    if not rejection and targetItem != -1:
        found = targetItem 
    else:
        found = -1
        
    return found
     
def determindMatch(inList, outList, myCountDic, folderName,contigReadGraph, N1):

    sd, newInList, newOutList = transformToCount(inList, outList, myCountDic, folderName,contigReadGraph, N1)

    return formResolveListForMatch(folderName, newInList, newOutList, sd, contigReadGraph, N1) 

def transformToCount(inList, outList, myCountDic, folderName,contigReadGraph, N1):

    newInList , newOutList = [], []
    
    newInList = getCt(inList ,myCountDic)
    newOutList = getCt(outList, myCountDic)
    
    sizeList = []
    for eachitem in myCountDic:
        sizeList.append(myCountDic[eachitem])
        
    sd = np.std(sizeList)
    
    return sd, newInList, newOutList

def transformToCountAggregate(inList, outList, myCountDic, folderName,contigReadGraph, N1,  Gnew, lenDic):

    newInList , newOutList = [], []
    
    newInList = getCtAgg(inList ,myCountDic, Gnew, lenDic)
    newOutList = getCtAgg(outList, myCountDic, Gnew, lenDic)
    
    sizeList = []
    for eachitem in myCountDic:
        sizeList.append(myCountDic[eachitem])
        
    sd = np.std(sizeList)
    
    return sd, newInList, newOutList

def formResolveListForMatch(folderName, newInList, newOutList, sd, contigReadGraph, N1):
    resolvedList = []

    for eachitem in newInList:
        found = satisfyMatch(eachitem ,newOutList, sd)
        
        if found != -1 :
            leftCtgIndex, rightCtgIndex = eachitem[0], found
            succReadsList = abunGraphLib.findPathBtwEnds(folderName, leftCtgIndex, rightCtgIndex, contigReadGraph, N1)
            
            if succReadsList != None:
                resolvedList.append([eachitem[0], found])
    
    return resolvedList

def determindMatchAggregate(inList, outList, myCountDic, folderName,contigReadGraph, N1, Gnew, lenDic):
    sd, newInList, newOutList = transformToCountAggregate(inList, outList, myCountDic, folderName,contigReadGraph, N1,  Gnew, lenDic)
    return formResolveListForMatch(folderName, newInList, newOutList, sd, contigReadGraph, N1) 

def addEdges(G, resolvedList):
    print "len(G.graphNodesList)",  len(G.graphNodesList)
    
    for eachitem in resolvedList:
        print eachitem
        G.insertEdge(eachitem[0], eachitem[-1], 1997)
    
def generateGapContentLookup(folderName, mummerLink, oldResolvedList, contigReadGraph, contigFilename,readsetFilename, mapDummyToRealDic={} ):
    gapContentLookUpList = []
    
    contigLenDic = IORobot.obtainLength(folderName, contigFilename+ ".fasta")
    N1 = len(contigLenDic)*2

    resolvedList  = []
    
    print "mapDummyToRealDic", mapDummyToRealDic
    for eachmatchpair in oldResolvedList:
        tmpList = []

        if eachmatchpair[0] >= N1 :
            tmpList = tmpList + mapDummyToRealDic[str(eachmatchpair[0] - N1)][1]
        else:
            tmpList.append(eachmatchpair[0])

        if eachmatchpair[-1] >= N1:
            tmpList = tmpList + mapDummyToRealDic[str(eachmatchpair[-1] - N1)][1]
        else:
            tmpList.append(eachmatchpair[-1])

        for ii in range(len(tmpList)-1):
            resolvedList.append([tmpList[ii], tmpList[ii+1]])
        
    gapContentLookUpList = abunGraphLib.parallelGapLookUp(resolvedList,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename)

    return gapContentLookUpList   

def abunSplit(folderName, mummerLink, myCountDic,contigReadGraph,  contigFilename,readsetFilename ):
    
    '''
    Input : repeatSpecification.txt , myCountDic.json, improved3.fasta, raw_reads.fasta
    Output : abunsplit.fasta
    
    Algorithm : 
    
    1. Load data from various sources [various json files]
    
    2. For each repeat interior:
        a) identify the abundances associated with in/out contigs
        b) perform a split and record the split
    
    3. Use split results to generate contigs [may already exist in newPhasing.py ] 
        a) use a graph to capture the split results 
        b) use reads to fill in any gaps 
        c) read out the contigs 
    
    '''
    
    
    json_data = open(folderName + "phaseRepeat.txt", 'r')
    repeatPairs = json.load(json_data)
    repeatPairs = obtainNonEmpty(repeatPairs)
    
    N1 = len(myCountDic)*2
    print "N1", N1
    
    G = graphLib.seqGraph(N1)
    
    gapContentLookUpList = []
    
    for eachitem in repeatPairs:
        inList, outList = eachitem[0], eachitem[1]
        resolvedList = determindMatch(inList, outList, myCountDic, folderName,contigReadGraph, N1)
        print "resolvedList", resolvedList
        gapContentLookUpList += generateGapContentLookup(folderName, mummerLink, resolvedList, contigReadGraph, contigFilename,readsetFilename)
        
        addEdges(G, resolvedList)
        
    gapContentLookUpDic = {}
    gapContentLookUpList.sort()
    
    
    for eachitem in gapContentLookUpList:
        gapContentLookUpDic[str(eachitem[0]) +"_" +str(eachitem[1])] = [ eachitem[2],eachitem[3],eachitem[4] ]
        print eachitem[2:4], len(eachitem[4])
    
    # some how change ASplitter here by appending necessary information   
    
    G.condense()
    IORobot.extractGraphToContigs(G, folderName, mummerLink, "abun.fasta", contigFilename+"_Double.fasta", gapContentLookUpDic)

def evaluateCoverage(dataList, lenDic, readLenDic, folderName,mummerLink, continueFilter, contigFilename):
    
    '''
    not sure if that is the right documentation... 

    Input : string_graph_3, improved3.fasta, raw_reads.fasta
    Output : string_graph_4 with weights [need a data structure to store the weight on node]

    Algorithm : 
    1. Find your favorite mappers to map read back
        a. MUMmer, Bowtie, bbmap, any that works V 
        b. And then write a short parser to parse the results V 
    2. Calculate count on the abundances 
        a. Aggregate by taking average [put weights on bin along contigs]
        b. Inheritance and a subclass 
    3. Find your favorite graphical tool to display 
        a. Use a javascript library [halfviz should just work ! put weight on edge ]

    '''    
    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = 0
    
    dataList.sort(key = itemgetter(-1)) 
    
    ctkk, ctbase = 0, 0
    toAddBackDic = copy.deepcopy(readLenDic)
    
    for key, items in groupby(dataList, itemgetter(-1)):
        maxMatch = -1
        bestname = ""
        
        for eachitem in items:
            ct = eachitem[6]/100.0 * eachitem[4]
            if ct > maxMatch:
                maxMatch = ct 
                bestname = eachitem[-2]
        myCountDic[bestname] += readLenDic[key] 
        
        ctkk = ctkk + 1 
        ctbase = ctbase + readLenDic[key]
        toAddBackDic[key] = -1
    
    cttot = 0
    for eachitem in readLenDic:
        cttot = cttot + readLenDic[eachitem]
        
    print "Missed coverage  ", (cttot - ctbase)/(4.7*pow(10, 6))
    print "percentage miss read", (len(readLenDic) - ctkk)/(1.0*len(readLenDic)) 
    
    toAddReadList = []
    for eachitem in toAddBackDic:
        if toAddBackDic[eachitem] >= 0 :
            toAddReadList.append(eachitem)
    
    '''
    This part need the most parallelism because it is most intense with -l 10 
    split V, workerList V , combine 
    '''
    
    if continueFilter:
        numberOfFiles= 20
        
        IORobot.putListToFileO(folderName, "raw_reads.fasta" , "selected_raw", toAddReadList)
        
        bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
        command = bindir + "/finisherSCCoreLib/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + "selected_raw.fasta"
        os.system(command)
        
        workerList = []
        
        for dummyI in range(1, numberOfFiles + 1):
            indexOfMum = ""
            if dummyI < 10:
                indexOfMum = "0" + str(dummyI)
            else:
                indexOfMum = str(dummyI)
           
            outputName, referenceName, queryName, specialName= "outAbunRefine"+indexOfMum, contigFilename+".fasta", "selected_raw.part-"+ indexOfMum + ".fasta",  "abunMissOut" + indexOfMum
            workerList.append([outputName, referenceName, queryName, specialName])
            
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,specialForRaw = True, refinedVersion = True)
        alignerRobot.combineMultipleCoorMum( True, mummerLink, folderName, "outAbunRefine", "abunMissOut", numberOfFiles)
        

        
    for eachitem in lenDic:
        #eachitem = "Segkk"+str(i)
        print eachitem , myCountDic[eachitem]/(1.0*lenDic[eachitem])
        myCountDic[eachitem] = myCountDic[eachitem]/(1.0*lenDic[eachitem])
        
    return myCountDic

def generateAbundanceGraph(folderName, mummerLink, contigFilename):
    
    
    print "generateAbundanceGraph"
    
    '''
    1. Find your favorite mappers to map read back
        a. MUMmer, Bowtie, bbmap, any that works V 
        b. And then write a short parser to parse the results V 
    '''
    numberOfFiles = 20
    workerList = []
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        
        '''
        "outGapFillRefine"+indexOfMum , "smaller_improvedContig.fasta",  "relatedReads_Double.part-" + indexOfMum + ".fasta",  "fromMumRefine" + indexOfMum
        '''
        outputName, referenceName, queryName, specialName= "outAbun"+indexOfMum, contigFilename+".fasta", "raw_reads.part-"+ indexOfMum + ".fasta",  "outAbun" + indexOfMum
        workerList.append([outputName, referenceName, queryName, specialName])
    
    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,False)
        '''
        command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + "out " + folderName + "improved3.fasta "+folderName+"raw_reads.part-" + indexOfMum + ".fasta"
        os.system(command)
    
        command = mummerLink + "show-coords -r " + folderName + "out.delta > " + folderName + "fromMumAbun" + indexOfMum
        os.system(command)
        '''
        
    dataList = []
    
    for i in range(1, 1+numberOfFiles): 
        if i < 10:
            indexOfMum = "0" + str(i)
        else:
            indexOfMum = str(i)
        dataList = dataList+ alignerRobot.extractMumData(folderName, "outAbun"+ str(indexOfMum)+"Out")
    

    '''
    2. Calculate count on the abundances 
        a. Aggregate by taking average [put weights on bin along contigs]
        b. Inheritance and a subclass 
    '''
         
    lenDic = IORobot.obtainLength(folderName, contigFilename+".fasta")
    readLenDic = IORobot.obtainLength(folderName , "raw_reads.fasta")
    

    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = [0 for i in range(lenDic[eachitem])]

    thres = 30
    lenSum = 0
    extraDataList= []
    
    
    print "len(dataList)", len(dataList)
    
    if not abunHouseKeeper.abunGlobalAvoidrefine: 
        myCountDic =  evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,  True, contigFilename)
        extraDataList = alignerRobot.extractMumData(folderName, "abunMissOut" )
    else:
        extraDataList = []
        
    dataList = dataList + extraDataList
    myCountDic = evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,False, contigFilename)
    
    with open(folderName + 'myCountDic.json', 'w') as f:
        json.dump(myCountDic, f)

    
    return myCountDic

def xNodeResolving(folderName, contigReadGraph):
    '''
    Input : contigGraph , abunInfo , folderName  

    Output: myresolvedList.json, gapContentLookUp.json, dummyNodeMapping.json

    Algorithm :
        1) Tranverse the graph 
            a) If the node can well be fixed with sd requirement met 
                i) Link it across and add the pair into the myresolvedList, gapContentLookUp
                ii) Add dummynodes and fill in the dummyNodeMapping 
        
        2) Format return and output as temp file 
    '''    
    
    ### Init G, myCountDic, N1
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)

    with open(folderName + 'myCountDic.json') as f:
        myCountDic = json.load(f)
    
    N1 = len(myCountDic)*2
    
    ### Add resolved edge  

    adj = [[] for i in range(N1)]
    
    for i in range(N1): 
        adj[i] = abunGraphLib.findAllReachable(i, N1, G)
    
    Gnew = graphLib.seqGraph(N1)
    
    for i in range(N1):
        for j in adj[i]:
            Gnew.insertEdge(i,j,1)
    
    extraCounter = 0
    mapDummyToRealDic = {}
    resolvedList = []
    
    for v in Gnew.graphNodesList:
        
        inList = [] 
        for eachitem in  v.listOfPrevNodes:
            inList.append(eachitem[0])
        
        outList = []
        for eachitem in v.listOfNextNodes:
            outList.append(eachitem[0])
        
        inListCt = getCtTwoToOne(inList,myCountDic)
        outListCt = getCtTwoToOne(outList,myCountDic)
    
        sizeList = []
        for eachitem in myCountDic:
            sizeList.append(myCountDic[eachitem])
            
        sd = np.std(sizeList)
        
        for eachIn in inListCt:
            matchedOut = satisfyMatch(eachIn , outListCt, sd) 
        
            if matchedOut != -1:
                leftCtgIndex, rightCtgIndex = eachIn[0], v.nodeIndex
                inSuccReadsList = abunGraphLib.findPathBtwEnds(folderName, leftCtgIndex, rightCtgIndex, contigReadGraph, N1)
                
                leftCtgIndex, rightCtgIndex =  v.nodeIndex, matchedOut
                outSuccReadsList = abunGraphLib.findPathBtwEnds(folderName, leftCtgIndex, rightCtgIndex, contigReadGraph, N1)
                
                if inSuccReadsList != None and outSuccReadsList != None :
                
                    resolvedList.append([eachIn[0] ]+inSuccReadsList+[ N1 + extraCounter])
                    print "in: ",  resolvedList[-1]
                    
                    resolvedList.append([N1+extraCounter]+outSuccReadsList+[ matchedOut])
                    print "out: ", resolvedList[-1]
                    
                    mapDummyToRealDic[extraCounter] = v.nodeIndex
                    extraCounter = extraCounter + 1

    return resolvedList, mapDummyToRealDic

def abunSplitWithXResolve(folderName, mummerLink, myCountDic,contigReadGraph,  contigFilename,readsetFilename ):
    N1 = len(myCountDic)*2
    print "N1", N1
    
    # Debug
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)

    adj = [[] for i in range(N1)]
    
    for i in range(N1): 
        adj[i] = abunGraphLib.findAllReachable(i, N1, G)
    
    Gnew = graphLib.seqGraph(N1)
    
    for i in range(N1):
        for j in adj[i]:
            Gnew.insertEdge(i,j,1)

    Gnew.reportEdge()
    # End Debug
    
    if False:
        json_data = open(folderName + "phaseRepeat.txt", 'r')
        repeatPairs = json.load(json_data)
        repeatPairs = obtainNonEmpty(repeatPairs)
            
        biResolvedCombineList = []
        for eachitem in repeatPairs:
            inList, outList = eachitem[0], eachitem[1]
            resolvedList = determindMatch(inList, outList, myCountDic, folderName,contigReadGraph, N1)
            
            biResolvedCombineList += resolvedList
        
        ### Xnode repeatResolution 
        xResolvedList, mapDummyToRealDic =  xNodeResolving(folderName, contigReadGraph)
        
        ### Combine resolution 
        resolvedList = xResolvedList + biResolvedCombineList
        resolvedList = abunHouseKeeper.getDistinct(resolvedList)
        print "resolvedList, len(resolvedList),len(xResolvedList), len(biResolvedCombineList) ", resolvedList, len(resolvedList), len(xResolvedList), len(biResolvedCombineList)
        
        
        with open(folderName + "resolvedList.json", 'w') as f:
            json.dump(resolvedList, f)
        
        with open(folderName + "mapDummyToRealDic.json", 'w') as f:
            json.dump(mapDummyToRealDic, f)
        

    if False:
        json_data = open(folderName + "resolvedList.json", 'r')
        resolvedList = json.load(json_data)
        
        json_data = open(folderName + "mapDummyToRealDic.json", 'r')
        mapDummyToRealDic = json.load(json_data)
        
        
        gapContentLookUpList = []
        gapContentLookUpList = generateGapContentLookup(folderName, mummerLink, resolvedList, contigReadGraph, contigFilename,readsetFilename, mapDummyToRealDic)
        gapContentLookUpDic = {}
        gapContentLookUpList.sort()
        
        for eachitem in gapContentLookUpList:
            gapContentLookUpDic[str(eachitem[0]) +"_" +str(eachitem[1])] = [ eachitem[2],eachitem[3],eachitem[4] ]
            print eachitem[2:4], len(eachitem[4])

        with open(folderName + "gapContentLookUpDic.json", 'w') as f:
            json.dump(gapContentLookUpDic, f)
        
    
    if False:
        json_data = open(folderName + "resolvedList.json", 'r')
        resolvedList = json.load(json_data)
        
        json_data = open(folderName + "mapDummyToRealDic.json", 'r')
        mapDummyToRealDic = json.load(json_data)
        
        G = graphLib.seqGraph(N1+ len(mapDummyToRealDic))
        addEdges(G, resolvedList)
        G.condense()
        
        G.saveToFile(folderName, "xResolvedGraph")
        


    if False:
        json_data = open(folderName + "mapDummyToRealDic.json", 'r')
        mapDummyToRealDic = json.load(json_data)
        
        G = graphLib.seqGraph(0)
        G.loadFromFile(folderName, "xResolvedGraph")
        
        json_data = open(folderName + "gapContentLookUpDic.json", 'r')
        gapContentLookUpDic = json.load(json_data)
            
        
        print "Final step: really hacking a file"
        os.system("cp "+ folderName + contigFilename+"_Double.fasta " +folderName + "tmpWithDummy.fasta")
        contigList = IORobot.readContigsFromFile(folderName,  contigFilename+"_Double.fasta")
        
        f = open(folderName + "tmpWithDummy.fasta", 'a')
        for i in range(len(mapDummyToRealDic)):
            id = mapDummyToRealDic[str(i)]
            f.write(">SegDum"+ str(i)+"\n")
            f.write(contigList[id] + "\n")
        f.close()
        
        
        IORobot.extractGraphToContigs(G, folderName, mummerLink, "abun.fasta", "tmpWithDummy.fasta", gapContentLookUpDic, mapDummyToRealDic)

def graphSurgery(myCountDic, folderName, contigReadGraph, mummerLink, readsetFilename,contigFilename ):

    ### Transitive reduction and remove double pointers
    N1 = len(myCountDic)*2
    print "N1", N1
    kthres = abunHouseKeeper.abunGlobalSplitParameterRobot.kthres
    edgeThres = abunHouseKeeper.abunGlobalSplitParameterRobot.edgeThres


    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, contigReadGraph)
    
    adj = [[] for i in range(N1)]
    
    for i in range(N1): 
        tmpList = abunGraphLib.findAllReachable(i, N1, G)
        
        for j in tmpList:
            if len(abunGraphLib.findAllPathK(i,j,G,kthres)) >= edgeThres:
                adj[i].append(j) 
    
    ### Filter adaptor skipped case 

    adaptorPair = []

    for i in range(len(adj)):
        if  i % 2 == 0:
            if i + 1 in adj[i]:
                adj[i].remove(i+1)
                adaptorPair.append([i, i+1])
        elif i % 2 ==1: 
            if i-1 in adj[i] :
                adj[i].remove(i-1)
                adaptorPair.append([i, i-1])
    
    Gnew = abunGraphLib.seqGraphDynamic(N1)
    
    for i in range(N1):
        for j in adj[i]:
            Gnew.insertEdge(i,j,1997)

    for eachpair in adaptorPair:
        u, v = eachpair[0], eachpair[1]
        for x in Gnew.graphNodesList[u].listOfPrevNodes:
            xIndex = x[0]
            Gnew.removeEdge(xIndex, v)
        for y in Gnew.graphNodesList[v].listOfNextNodes:
            yIndex = y[0]
            Gnew.removeEdge(u, yIndex)


    Gnew.reportEdge()
    ### Trying out the new component 
    import toCondenseFixer
    Gnew = toCondenseFixer.noGoZoneDefiner(Gnew, folderName)

    Gnew.symGraph()
    #Gnew.reportEdge()
    ### End filter adaptor skipped case 

    if abunHouseKeeper.abunGlobalSplitParameterRobot.runGraphSurgery:
    
        Gnew.initAdv()    
        if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunCondenseRemove:
            Gnew.condenseEdgeRemove(G, folderName, mummerLink, contigFilename)

        if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunDoubltPtr:
            Gnew.doubleEdgeReduction()
        
        if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunTransitive:
            Gnew.transitiveReduction(folderName, mummerLink, contigFilename+"_Double.fasta", readsetFilename+"_Double.fasta", G)
        
        Gnew.condense()
        Gnew.findAdjList()
    else:
        Gnew.initAdv()    
        Gnew.condense()
        Gnew.findAdjList()
    
    return Gnew

def BResolution(Gnew, folderName, contigReadGraph, N1, myCountDic, lenDic, mummerLink):

    if abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve:
        print "abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve", abunHouseKeeper.abunGlobalSplitParameterRobot.runBResolve
        maxRThres = abunHouseKeeper.abunGlobalSplitParameterRobot.RThres

        repeatFinder.adjListToRepeatList(Gnew.adj,folderName,"phaseRepeatTR.txt")

        json_data = open(folderName + "phaseRepeatTR.txt", 'r')
        repeatPairs = json.load(json_data)


        repeatPairs = obtainNonEmpty(repeatPairs)
        
        biResolvedCombineList = []

        G = abunGraphLib.seqGraphWt(0)
        G.loadFromFile(folderName, contigReadGraph)

        Grev = abunGraphLib.formReverseGraphFast(G)

        abunAnalysisList = []

        for eachitem in repeatPairs:
            inList, outList = eachitem[0], eachitem[1]
            if not abunHouseKeeper.abunGlobalRunEM:
                resolvedList, brResolvedList = [], [] 

                if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunAbunB:
                    if abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerB > 0:
                        abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLower = abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerB

                    if abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperB > 0:
                        abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpper = abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperB

                    if not abunHouseKeeper.abunGlobalSplitParameterRobot.toRunAggB: 
                        resolvedList = determindMatch(inList, outList, myCountDic, folderName,contigReadGraph, N1)
                    else:

                        resolvedList = determindMatchAggregate(inList, outList, myCountDic, folderName,contigReadGraph, N1, Gnew, lenDic)
                                    

                if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunBRB:
                    if abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresB > 0:
                        abunHouseKeeper.abunGlobalSplitParameterRobot.BRThres = abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresB

                    brResolvedList = formBRReolve(folderName, inList, outList, G, Grev, True, N1)


                combinedList = abunHouseKeeper.getDistinct(resolvedList + brResolvedList)
                
                print "resolvedList, brResolvedList, inList, outList", resolvedList, brResolvedList, inList, outList
                
                print "resolveConflict(combinedList)", resolveConflict(combinedList)      

                abunAnalysisList.append([inList, outList,resolvedList, brResolvedList, resolveConflict(combinedList) ])
                if  len(inList) <= maxRThres and  len(outList) <= maxRThres and len(inList) > 0 and len(outList) > 0:
                    resolvedCombine = resolveConflict(combinedList)
                    Gnew.bipartiteLocalResolve(resolvedCombine , inList, outList, folderName)
            else:
                resolvedCombine = emalgo.BResolvePreparation(folderName, inList, outList,  G, Grev, N1, mummerLink)
                Gnew.bipartiteLocalResolve(resolvedCombine , inList, outList, folderName)

        Gnew.condense()

        with open(folderName + "biResolvedCombineList.json", 'w') as f:
            json.dump(biResolvedCombineList, f)    


        with open(folderName + "abunAnalysisList.json", 'w') as f:
            json.dump(abunAnalysisList, f)    

        #assert(1==2)


        return Gnew

    else:
        return Gnew
            
def XResolution(folderName,contigReadGraph, Gnew, myCountDic, lenDic, N1, mummerLink):

    if abunHouseKeeper.abunGlobalSplitParameterRobot.runXResolve:
        G = graphLib.seqGraph(0)
        G.loadFromFile(folderName, contigReadGraph)

        Grev = abunGraphLib.formReverseGraphFast(G)
        
        
        
        if not abunHouseKeeper.abunGlobalRunEM:
            xResolvedList, brResolvedListforX = [[] for i in range(N1)], [[] for i in range(N1)]
            if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunAbunX:
                if abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerX > 0:
                    abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLower = abunHouseKeeper.abunGlobalSplitParameterRobot.AbunLowerX

                if abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperX > 0:
                    abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpper = abunHouseKeeper.abunGlobalSplitParameterRobot.AbunUpperX

                xResolvedList =  xNodeAdvResolving(Gnew, G, folderName, myCountDic, lenDic)
                
            if abunHouseKeeper.abunGlobalSplitParameterRobot.toRunBRX:
                
                if abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresX > 0:
                    abunHouseKeeper.abunGlobalSplitParameterRobot.BRThres = abunHouseKeeper.abunGlobalSplitParameterRobot.BRThresX

                brResolvedListforX = xNodeBrResolving(Gnew, G,Grev, folderName, N1)

            combinedList = resolveConflictX(xResolvedList, brResolvedListforX) 
        else:
            combinedList = xNodeEMResolving(Gnew, G, Grev, folderName, myCountDic, lenDic, N1, mummerLink)
        
        print "combinedList", combinedList

        Gnew.xResolve(combinedList)        
        Gnew.condense()
        Gnew.saveToFile(folderName, "xResolvedGraph")
        
        with open(folderName + "mapDummyToRealDic.json", 'w') as f:
            json.dump(Gnew.mapDummyToRealDic, f)    
        
        with open(folderName + "xResolvedSimplifiedList.json", 'w') as f:
            json.dump(Gnew.xResolvedSimplifiedList, f)    


    else:
        Gnew.saveToFile(folderName, "xResolvedGraph")
        
        with open(folderName + "mapDummyToRealDic.json", 'w') as f:
            json.dump(Gnew.mapDummyToRealDic, f)    
        
        with open(folderName + "xResolvedSimplifiedList.json", 'w') as f:
            json.dump(Gnew.xResolvedSimplifiedList, f)   

def formGapLookUp(folderName,contigReadGraph, mummerLink,contigFilename,readsetFilename):

    json_data = open(folderName + "biResolvedCombineList.json", 'r')
    biResolvedCombineList = json.load(json_data)
    
    json_data = open(folderName + "xResolvedSimplifiedList.json", 'r')
    xResolvedSimplifiedList = json.load(json_data)
    
    json_data = open(folderName + "mapDummyToRealDic.json", 'r')
    mapDummyToRealDic = json.load(json_data)
    
    resolvedList = biResolvedCombineList + xResolvedSimplifiedList
    
    gapContentLookUpList = []
    gapContentLookUpList = generateGapContentLookup(folderName, mummerLink, resolvedList, contigReadGraph, contigFilename,readsetFilename, mapDummyToRealDic)
    gapContentLookUpDic = {}
    gapContentLookUpList.sort()
    
    for eachitem in gapContentLookUpList:
        gapContentLookUpDic[str(eachitem[0]) +"_" +str(eachitem[1])] = [ eachitem[2],eachitem[3],eachitem[4] ]
        print eachitem[2:4], len(eachitem[4])

    with open(folderName + "gapContentLookUpDic.json", 'w') as f:
        json.dump(gapContentLookUpDic, f)

def readContigForAbunSplit(folderName,mummerLink,  contigFilename, readsetFilename, N1, contigReadGraph):

    json_data = open(folderName + "mapDummyToRealDic.json", 'r')
    mapDummyToRealDic = json.load(json_data)
    
    G =[] 
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, "xResolvedGraph")
    
    gapContentLookUpDic = {}

    furtherGapList = []
    for i in range(N1):
        if len(G.graphNodesList[i].nodeIndexList) > 1:
            for j in range(len(G.graphNodesList[i].nodeIndexList)-1):
                
                bk, fwd = G.graphNodesList[i].nodeIndexList[j] , G.graphNodesList[i].nodeIndexList[j+1]
                
                key = str(bk)+ "_"+ str(fwd)
                
                if not key in gapContentLookUpDic:
                    furtherGapList.append([bk,fwd])
     
    with open(folderName + "furtherGapList.json", 'w') as f:
        json.dump(furtherGapList, f)
        
    
    furtherGapContentLookUpList = generateGapContentLookup(folderName, mummerLink, furtherGapList , contigReadGraph, contigFilename,readsetFilename, mapDummyToRealDic)
    
    
    for eachitem in furtherGapContentLookUpList:
        gapContentLookUpDic[str(eachitem[0]) +"_" +str(eachitem[1])] = [ eachitem[2],eachitem[3],eachitem[4] ]
        print eachitem[2:4], len(eachitem[4])
    
    #segLookUp = IORobot.readContigsFromFile(folderName, "LC_n_Double.fasta")
    
    print "Final step: really hacking a file"
    os.system("cp "+ folderName + contigFilename+"_Double.fasta " +folderName + "tmpWithDummy.fasta")
    contigList = IORobot.readContigsFromFile(folderName,  contigFilename+"_Double.fasta")
    
    IORobot.extractGraphToContigs(G, folderName, mummerLink, "abunPre.fasta", "tmpWithDummy.fasta", gapContentLookUpDic, mapDummyToRealDic)

    if True:
        nonRedundantResolver.removeRedundantWithFile(folderName , mummerLink, "abunPre", "abunMum", "abun")

def xNodeEMResolving(Gnew, GContigRead,Grev, folderName, myCountDic, lenDic, N1, mummerLink):
    print "emalgo"
    combinedList = emalgo.XResolvePreparation(Gnew, GContigRead,Grev, folderName, myCountDic, lenDic, N1, mummerLink)
    print "combinedList", combinedList
    return combinedList
    
def xNodeAdvResolving(Gnew, GContigRead, folderName, myCountDic, lenDic):
    N1 = len(myCountDic)*2
    
    xResolvedList = [ [] for i in range(N1)]
    
    for nodeI in range(N1):
        
        v = Gnew.graphNodesList[nodeI]
        
        inList = [] 
        for eachitem in  v.listOfPrevNodes:
            inList.append(eachitem[0])
        
        outList = []
        for eachitem in v.listOfNextNodes:
            outList.append(eachitem[0])
        

        if not abunHouseKeeper.abunGlobalSplitParameterRobot.toRunAggX:
            inListCt  =   getCtTwoToOne(inList ,myCountDic)
            outListCt =   getCtTwoToOne(outList,myCountDic)
        else:
            inListCt  =   getCtTwoToOneAgg(inList ,myCountDic, Gnew, lenDic)
            outListCt =   getCtTwoToOneAgg(outList,myCountDic, Gnew, lenDic)
    
        sizeList = []
        for eachitem in myCountDic:
            sizeList.append(myCountDic[eachitem])
            
        sd = np.std(sizeList)

        print "SD: " , sd, len(inListCt), len(inList)
        
        for eachIn in inListCt:
            matchedOut = satisfyMatch(eachIn , outListCt, sd) 
            print "matchedOut", matchedOut

            if matchedOut != -1:
                leftCtgIndex, rightCtgIndex = eachIn[0], v.nodeIndex
                inSuccReadsList = abunGraphLib.findPathBtwEndsFast(folderName, leftCtgIndex, rightCtgIndex, GContigRead, N1)
                
                leftCtgIndex, rightCtgIndex =  v.nodeIndex, matchedOut
                outSuccReadsList = abunGraphLib.findPathBtwEndsFast(folderName, leftCtgIndex, rightCtgIndex, GContigRead, N1)
                
                if inSuccReadsList != None and outSuccReadsList != None :                
                    xResolvedList[nodeI].append([eachIn[0], matchedOut ])
                    
    return xResolvedList

def formBRReolve(folderName, oldin, oldout, G, Grev, needTransform, N1):
    inList, outList = [],[]
    
    if needTransform : 
        for eachitem in oldin : 
            inList.append(eachitem/2)
        for eachitem in oldout : 
            outList.append(eachitem/2)
    else:
        for eachitem in oldin : 
            inList.append(eachitem)
        for eachitem in oldout : 
            outList.append(eachitem)

    resolvedList = formConfirmReadResolve(folderName, inList, outList, G, Grev,  N1)

    newResolvedList = filterConfidResolve(resolvedList)

    print "newResolvedList", newResolvedList 
    return newResolvedList 

def filterConfidResolve(resolvedList):
    newResolvedList = []
    resolvedList.sort()
    conThres = abunHouseKeeper.abunGlobalSplitParameterRobot.BRThres
    print "conThres", conThres

    for key, items in groupby(resolvedList):
        tmpList = list(items)
        if len(tmpList) >= conThres:
            newResolvedList.append(key)

    if False:
        noConflict = resolveConflict(abunHouseKeeper.getDistinct(resolvedList))
        noConflict = abunHouseKeeper.getDistinct(noConflict)
        newResolvedList = abunHouseKeeper.getDistinct(newResolvedList)
        newResolvedList = abunHouseKeeper.getDistinct(intersect(newResolvedList , noConflict))

    return newResolvedList

def intersect(a, b):
    c = a +b 
    c.sort()
    returnList = []
    for i in range(len(c)-1):
        if c[i] == c[i+1]:
            returnList.append(c[i])
            
    return returnList 

def resolveConflict(combinedList):
    combinedDic = {}
    
    for eachitem in combinedList:
        combinedDic[str(eachitem[0])+"_"+str(eachitem[1])] = True

    for i in range(2):
        combinedList.sort(key = itemgetter(i))

        for key, items in groupby(combinedList, itemgetter(i)):
            tmpList = list(items)
            if len(tmpList) > 1: 
                for eachitem in tmpList:
                    combinedDic[str(eachitem[0])+"_"+str(eachitem[1])] = False

    resolvedList = []

    combinedList.sort()
    for key, items in groupby(combinedList):
        name = str(key[0]) + "_" + str(key[1])
        if combinedDic[name] == True: 
            resolvedList.append(key)

    return resolvedList 
 
def formConfirmReadResolve(folderName, inList, outList, G, Grev,  N1):
    #print "formConfirmReadResolve"

    resolvedList = []
    confirmingReadList = [] 
    brLFlankList = []
    brRFlankList = []


    ### Find possible candidate reads 
    print "inList , outList formConfirmReadResolve()", inList, outList 
    for eachin in inList:
        for eachout in outList:
            pathList = abunGraphLib.findAllPathK(eachin,eachout,G,3)
            for path in pathList:
                if len(path) == 3 and path[1]>= N1 : 
                    R = path[1] 
                    confirmingReadList.append(R)
                    brLFlankList.append([eachin, R])
                    brRFlankList.append([eachout, R])


    ### Filter simple false cases
    toUseReadDic = {}
    confirmingReadList.sort()
    for key, items in groupby(confirmingReadList):
        toUseReadDic[str(key)] = True


    newbrLFlankList = abunHouseKeeper.getDistinct(brLFlankList)
    newbrLFlankList.sort(key = itemgetter(1))

    for key, items in groupby(newbrLFlankList, itemgetter(1)):
        mylist = list(items)
        if len(mylist) > 1: 
            toUseReadDic[str(key)] = False

    newbrRFlankList = abunHouseKeeper.getDistinct(brRFlankList)
    newbrRFlankList.sort(key = itemgetter(1))
    
    for key, items in groupby(newbrRFlankList, itemgetter(1)):
        mylist = list(items)
        if len(mylist) > 1: 
            toUseReadDic[str(key)] = False

    finalSearchReadList = []
    for eachitem in toUseReadDic:
        if toUseReadDic[eachitem] == True:
            finalSearchReadList.append(int(eachitem))


    ### Check paths to confirm all false cases
    for eachR in finalSearchReadList:
        l1 = abunGraphLib.findAllReachable(eachR, N1, G)
        l2 = abunGraphLib.findAllReachable(eachR, N1, Grev)

        l1Distinct = abunHouseKeeper.getDistinct(l1)
        l2Distinct = abunHouseKeeper.getDistinct(l2)

        if len(l1Distinct) == 1 and len(l2Distinct) ==1 : 
            c1, c2 = l1Distinct[0], l2Distinct[0]
            resolvedList.append([c2, c1])

    return resolvedList

def resolveConflictX(listA, listB):
    resolvedList = [[] for i in range(len(listA))]
    print "len(listA), len(listB)", len(listA), len(listB)
    for i in range(len(listA)):
        combinedList = listA[i] + listB[i]
        newCombinedList = abunHouseKeeper.getDistinct(combinedList)
        tmpResolved = resolveConflict(newCombinedList)
        resolvedList[i] = tmpResolved 

    return resolvedList

def xNodeBrResolving(Gnew, G, Grev, folderName, N1):
    #brResolvedListforX = []
    xResolvedList = [ [] for i in range(N1)]
    
    for nodeI in range(N1):
        
        v = Gnew.graphNodesList[nodeI]
        
        inList = [] 
        for eachitem in  v.listOfPrevNodes:
            inList.append(eachitem[0])
        
        outList = []
        for eachitem in v.listOfNextNodes:
            outList.append(eachitem[0])

        xResolvedList[nodeI] = formBRReolve(folderName, inList, outList, G, Grev, False, N1)

    return xResolvedList 

def abunSplitAdvResolve(folderName, mummerLink, myCountDic,contigReadGraph,  contigFilename,readsetFilename ):

    '''
    Algorithm: 
    1)Load ContigReadGraph and form xResolvedGraph
    2)Transitive reduction and remove double pointers
    3)Bipartite resolution
    4)xResolve 
    5)Form gapLookUp 
    6)Read contigs out from graph
    7)CheckAns and get it done today again... 
    '''
    emalgo.generateAssociatedReadDic(folderName) 
    
    lenDic = IORobot.obtainLength(folderName, contigFilename+"_Double.fasta")
    N1 = len(lenDic)

    Gnew = graphSurgery(myCountDic, folderName, contigReadGraph, mummerLink, readsetFilename, contigFilename)
    Gnew.logEdges(folderName, "graphsurgery")
    
    #Gnew.reportEdge()
    #assert(False)

    Gnew = BResolution(Gnew, folderName, contigReadGraph, N1, myCountDic, lenDic, mummerLink)
    Gnew.logEdges(folderName, "BResolution")
    
    XResolution(folderName,contigReadGraph, Gnew, myCountDic, lenDic, N1 , mummerLink)
    Gnew.logEdges(folderName, "XResolution")

    readContigForAbunSplit(folderName,mummerLink,  contigFilename, readsetFilename, N1,contigReadGraph)

def splitter(folderName, mummerLink, contigReadGraph, contigFilename,readsetFilename ):

    '''
    Input : repeatSpecification.txt , myCountDic.json, improved3.fasta, raw_reads.fasta
    Output : abunsplit.fasta
    
    Algorithm : 
    
    1. Load data from various sources [various json files]
    
    2. For each repeat interior:
        a) identify the abundances associated with in/out contigs
        b) perform a split and record the split
    
    3. Use split results to generate contigs [may already exist in newPhasing.py ] 
        a) use a graph to capture the split results 
        b) use reads to fill in any gaps 
        c) read out the contigs 
    
    '''
    
    
    with open(folderName + 'myCountDic.json') as f:
        myCountDic = json.load(f)
        
    #abunSplit(folderName, mummerLink, myCountDic, contigReadGraph, contigFilename,readsetFilename )    
    #abunSplitWithXResolve(folderName, mummerLink, myCountDic,contigReadGraph,  contigFilename,readsetFilename )
    abunSplitAdvResolve(folderName, mummerLink, myCountDic,contigReadGraph,  contigFilename,readsetFilename )
    
def mainFlow(folderName, mummerLink):
    print "Hello world"
    
    contigFilename = "improved3"
    readsetFilename = "phasingSeedName"
    optTypeFileHeader = "phaseString"
    contigReadGraph = "phaseStringGraph1"
    repeatFilename = "phaseRepeat.txt"
    repeatSpec = "repeatSpecification.txt"
    optionToRun = "xphase"
    
        
    if abunHouseKeeper.abunGlobalRunPickUp == "map" :
        associatedReadFinder.getAllAssociatedReads(folderName, mummerLink,readsetFilename)
        readContigGraphFormer.formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph )
        repeatFinder.identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun )
        
    
    if abunHouseKeeper.abunGlobalRunPickUp == "map" or abunHouseKeeper.abunGlobalRunPickUp == "count" :
        myCountDic = generateAbundanceGraph(folderName, mummerLink, contigFilename)
        
    if abunHouseKeeper.abunGlobalRunPickUp == "map" or abunHouseKeeper.abunGlobalRunPickUp == "count" or abunHouseKeeper.abunGlobalRunPickUp == "split" :
        splitter(folderName, mummerLink, contigReadGraph, contigFilename,readsetFilename )

    if abunHouseKeeper.abunGlobalRunPickUp == "graph":
        print "Graph here"
        readContigGraphFormer.formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph, False)
        splitter(folderName, mummerLink, contigReadGraph, contigFilename,readsetFilename )
    
        
    os.system("cp selected_raw.part-* "+ folderName )
    os.system("rm selected_raw.part-*")
  















