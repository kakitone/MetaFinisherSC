from ..repeatPhaserLib.finisherSCCoreLib import IORobot
from ..repeatPhaserLib.finisherSCCoreLib import alignerRobot
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver 
import bisect


mergerGlobalLCReads = "SR"


import os 
import random 
import json 
from operator import itemgetter
from itertools import groupby 
import numpy as np


'''
Goal : Merging contigs from two libraries in a robust way 

High level interface  : 
    Input: SR.fasta, LR.fasta, SC.fasta, LC.fasta 
    Output: contigs.fasta
    Algorithm: 
        1. Obtain a mis-assembly fixed and non-redundant contigs from LC; 
            call it LC_n
        2. Obtain a mis-assembly fixed and non-redundant contigs from SC; 
            call it SC_n
        3. Combine SC_n and LC_n and output it as contigs.fasta
'''
    
    
def combineResults(folderName):
    command = "cat "+ folderName + "LC_n.fasta "+ folderName + "SC_n.fasta > "
    command = command + folderName + "contigs.fasta"  
    
    if True:  
        os.system(command)
        
def mergeContigs(folderName, mummerLink):
    print "mergeContigs" 
    
    fixLCMisassembly(folderName, mummerLink)
    fixSCMisassembly(folderName, mummerLink)
    
    combineResults(folderName)
    
'''        

Subroutine 1: Obtain LC_n
    Input : SR.fasta, LR.fasta, LC.fasta
    Output : LC_n.fasta
    Algorithm : 
        1. Align reads to contigs
        2. Calculate abundance information and make break at 
            suspicious junctions
    

    Remark : criterion to break == imbalance of abundances

'''




def completelyEmbedSR2TR(eachsub, lenDic) :
    #  "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    
    if mergerGlobalLCReads == "SR":        
        thres = 7
        
        readName = eachsub[-1]
        matchReadLen = eachsub[5]
        
        if abs(lenDic[readName] - matchReadLen ) < thres:
            return True
        else:
            return False
    elif mergerGlobalLCReads == "LR":
        thres = 20
        rdStart , rdEnd = min(eachsub[2], eachsub[3]), max(eachsub[2], eachsub[3])
        
        readName = eachsub[-1]
        readLen = lenDic[readName] 
        
        if rdStart < thres and rdEnd >= readLen - thres:
            return True
        else:
            return False
        

def assignCoverage(dataitem, coveragePerContigsDic):
    contigName = dataitem[-2]
    startPt, endPt = dataitem[0], dataitem[1]
    
    if not contigName in coveragePerContigsDic:
        print "contig not found" + contigName
    
    else:
        if endPt > len(coveragePerContigsDic[contigName]) :
            print "End point too far", endPt
        if startPt < 0 :
            print "start point too far", startPt
            
    for i in range(startPt-1, endPt):
        coveragePerContigsDic[contigName][i] += 1
    
    

def alignSR2LC(folderName, mummerLink):
    print "alignSR2LC"
    '''
    Input : SR.fasta, LC.fasta 
    Output : 
        coveragePerContigs.json : [[50, 23, 41] , [11, 23 ,45,77] , ...]
        Aka coverage plot for each contigs 
    
    Algorithm : 
        V 0) Seed random with the same number to insist consistency
        V 1) Use MUMmer align to map short reads to long contigs. If multiple matches are detected with certain tolerance, then randomly assign to one of the copies. 
        V 2) Evaluate count 
        V 3) Save to file  coveragePerContigs.json
        
    Rmk: Beaware of double stranded DNA
        
    '''
    LCLenDic = IORobot.obtainLength(folderName, "LC.fasta")
    coveragePerContigsDic = {}
    for eachitem in LCLenDic : 
        coveragePerContigsDic[eachitem] = [ 0 for i in range(LCLenDic[eachitem]) ]
    
    random.seed(0)
    outputName, referenceName, queryName = "SR2LCAlign" , "LC", "SR" 
    
    if True:
        dataList = alignerRobot.largeRvsQAlign(folderName, 20, mummerLink, referenceName, queryName, outputName)

    
    dataList.sort(key = itemgetter(-1))
    LCLenDic = IORobot.obtainLength(folderName, "LC.fasta")
    SRLenDic = IORobot.obtainLength(folderName, "SR.fasta")
    
    
    for eachitem in coveragePerContigsDic:
        print eachitem
        
    for key, items in groupby(dataList, itemgetter(-1)) :
        itemList = [] 
        for eachsub in items:
            if completelyEmbedSR2TR(eachsub, SRLenDic) :
                itemList.append(eachsub)
        
        if len(itemList) > 0:
            i = random.randint(0,len(itemList)-1)
            assignCoverage(itemList[i], coveragePerContigsDic)
        
    with open(folderName + "coveragePerContigs.json", 'w') as outfile:
        json.dump(coveragePerContigsDic, outfile)
    

def obtainTDiffVec(covVec, T):
    n = len(covVec)
    tDiffVec = [0 for i in range(n)]
    
    
    lastBack = np.mean(covVec[0:T])
    lastFor = np.mean(covVec[T+1:T+T+1])
    
    tDiffVec[T] = lastFor- lastBack
    
    for i in range(T+1, n-T):
        #bkwd = np.mean(covVec[i-T:i])
        #fwd = np.mean(covVec[i+1:i+T+1])
        bkwd = lastBack + 1.0/T*(covVec[i-1] - covVec[i-T-1])
        fwd = lastFor + 1.0/T*(-covVec[i] + covVec[i+T] )
        
        delta = fwd - bkwd
        tDiffVec[i] = delta
        
        lastBack = bkwd
        lastFor = fwd
        
    return tDiffVec 


def interQuartile_sd(tDiffVec, T):
    sd = 1
    quartile = np.percentile(tDiffVec, 25)
    thirdQuartile = np.percentile(tDiffVec, 75)
    stdList = []
    
    n = len(tDiffVec)
    
    for i in range(T, n-T):
        if quartile <= tDiffVec[i] <= thirdQuartile:
            stdList.append(tDiffVec[i])  
    
    sd = np.std(stdList)
    
    return sd

    
def findOutliners(tDiffVec, sd, T):
    outLiners = []
    
    n = len(tDiffVec)
    deletedList = [False for i in range(n)]
    searchList = [[abs(tDiffVec[i]), i ] for i in range(n)]
    searchList.sort(reverse = True)
    print searchList[0:10]
    
    for eachitem in searchList:
        value, index = eachitem[0] , eachitem[1]
        #print value, index
        if value < 13*sd :
            break
        else:
            if deletedList[index] == True:
                continue
            else:
                outLiners.append(index)
                for j in range(index - T,  index+ T+1):
                    deletedList[j] = True
    
    return outLiners


def filterBreakPts(outLiners, T, covVec, sd):
    modifiedOutliners =  []
    x = sorted(outLiners)
    m = len(outLiners)
    n = len(covVec)
    
    x.insert(0, 0)
    x.insert(m+1, n-1)
    
    Tnew = T/4
    
    for i in range(1, m+1):
        begin, end = x[i-1] + Tnew , x[i] - Tnew
        bkwd= np.mean(covVec[begin:end])
        
        begin, end = x[i] + Tnew, x[i+1] - Tnew
        fwd = np.mean(covVec[begin:end])
        
        if abs(bkwd - fwd) > 5*sd:
            modifiedOutliners.append(x[i])
    
    m = len(modifiedOutliners)
    
    modifiedOutliners.insert(0, 0)
    modifiedOutliners.insert(m+1, n-1)
            
    return modifiedOutliners  

def fineToneBreakPts(outLiners, T, covVec):
### Deprecated... not really needed
    modifiedOutliners =  []
    x = sorted(outLiners)
    m = len(outLiners)
    n = len(covVec)
    
    x.insert(0, 0)
    x.insert(m+1, n-1)
    
    Tnew = T/4
    
    for i in range(1, m+1):
        begin, end = x[i-1] + Tnew , x[i] - Tnew
        bkwd= np.mean(covVec[begin:end])
        
        begin, end = x[i] + Tnew, x[i+1] - Tnew
        fwd = np.mean(covVec[begin:end])
        
        obj = 0
        begin , end = x[i] - Tnew , x[i] + Tnew 
        for t in range(begin, end +1):
            obj = obj+ ( covVec[t] -fwd)**2
        
        minVal , index = obj, x[i] - Tnew -1
        
        for t in range(begin, end):
            delta =  (covVec[t] - bkwd)**2  - (covVec[t] - fwd)**2
            obj = obj + delta
            print delta, bkwd, fwd, covVec[t]
            if obj < minVal:
                minVal = object
                index = t
                
        modifiedOutliners.append(index)
        print "index, begin, end", index, begin, end
    
    return modifiedOutliners


def breakAcBkPts(contig, modifiedOutliners):
    contigBreakDown = []
    m = len(modifiedOutliners) - 2
    x = sorted(modifiedOutliners)
    
    for i in range(m+1):
        start , end = x[i], x[i+1]
        newItem = contig[start:end]
        contigBreakDown.append(newItem)
        
    return contigBreakDown


    

def breakAcBkPtsTwoSided(contig, modifiedOutliners, folderName, mummerLink):
    contigBreakDown = []
    
    modifiedOutliners.sort()
    
    segList = []
    
    for i in range(len(modifiedOutliners)-1):
        
        startPt = modifiedOutliners[i] 
        endPt = modifiedOutliners[i+1]
        
        # New version : breakDic format : 
        # {'seg1' : [ [contigName, bkPtL, bkPtR, readNameL, readNameR] ... ], [[],[], ... ]}
        # insertItem = [begin, end, leftRead, rightRead]

        insertItem = []
        
        if startPt[1] != None and endPt[2] != None:
            insertItem = [startPt[1], endPt[2], startPt[3], endPt[4] ]
        elif startPt[1] == None :
            insertItem = [startPt[2], endPt[2], None, endPt[4] ]
        elif endPt[2] == None :
            insertItem = [startPt[1], endPt[1], startPt[3], None ]
             
        segList.append(insertItem)
        
    for eachitem in segList:
        leftpart, middlepart, rightpart = "", "", ""
        begin, end, leftRead, rightRead = eachitem[0],eachitem[1],eachitem[2],eachitem[3]
        print "begin, end, leftRead, rightRead: ", begin, end, leftRead, rightRead
        if leftRead != None:
            leftSeg, rightSeg = IORobot.myRead(folderName, "LR.fasta",leftRead ), contig[begin:end] 
            overlap = IORobot.align(leftSeg, rightSeg, folderName, mummerLink)
            print "overlap:", overlap
            
            leftpart = leftSeg[0:len(leftSeg)-overlap[0]]
          
        if rightRead != None:
            leftSeg, rightSeg = contig[begin:end], IORobot.myRead(folderName, "LR.fasta",rightRead ) 
            overlap = IORobot.align(leftSeg, rightSeg, folderName, mummerLink)
            print "overlap:", overlap
            
            rightpart = rightSeg[overlap[1]:]
        
        #assert(1==2)
        middlepart = contig[begin:end] 
        
        combinedSeg = leftpart + middlepart + rightpart
        print "len(combinedSeg)", len(combinedSeg)
        contigBreakDown.append(combinedSeg)
        
    return contigBreakDown 

def breakLC(folderName):
    
    '''
    Input : LC.fasta, coveragePerContigs.json
    Output : LC_n.fasta
    Algorithm : 
        V 0) Set T  = 1K ; contigList = [] ; Read in LC and coveragePerContigs
        For each contig x in LC.fasta: 
            1) Obtain first difference of length T vector y
            2) Calculate sd = mod_sd(x)
            3) Find points with value > 3sd by
                a) sort
                repeat 
                    b) take the highest
                    c) if > 3sd continue ; else break
                    d) disable +/-T values
                until no satisfying items
            4) Point point the breaking point by using longer averages and l2 min 
            5) If there is any breaking locations, break x into pieces and save it broken list to contigList
        V 6) Write contigList to LC_n.fasta
    
    Remark: Think a bit more ... looks like a batch to be deleted ? !
    '''
    T = 1000 
    multThres = 12
    
    contigList = []
    
    json_data = open(folderName + "coveragePerContigs.json", 'r')
    coveragePerContigs = json.load(json_data)
    
    LCList = IORobot.readContigsFromFile(folderName, "LC.fasta")
    lenDic = IORobot.obtainLength(folderName, "LC.fasta")
    
    name2Index, index2Name = IORobot.fastaContigNameIndexConversion(folderName, "LC.fasta")
    
    for contigIndex in range(len(LCList)) : 
        
        name = index2Name[contigIndex]
        covVec = coveragePerContigs[name]
        
        if lenDic[name] > multThres*T:
            tDiffVec = obtainTDiffVec(covVec, T)
            print "len(tDiffVec)", len(tDiffVec)
            
            sd = interQuartile_sd(tDiffVec, T)
            print "sd", sd
            
            outLiners = findOutliners(tDiffVec, sd, multThres/2*T)
            print "outLiners", outLiners
            
            modifiedOutliners = filterBreakPts(outLiners, T, covVec, sd)
            print "modifiedOutliners", modifiedOutliners
            
            # Note : seems like not needed, so not use
            # modifiedOutliners = fineToneBreakPts(outLiners, T, covVec)
            
            contigBreakDown = breakAcBkPts(LCList[contigIndex], modifiedOutliners)
            
            contigList = contigList + contigBreakDown
            #assert(1==2)
        else:
            contigList = contigList + [LCList[contigIndex]]
            
    
    IORobot.writeSegOut(contigList, folderName, "LC_n.fasta")
    

def fixLCMisassembly(folderName, mummerLink):
    print "fixLCMisassembly"
    # Remark : if done, one can also try on LR only data set to see if there 
    # is any reduction on mis-assembly. 
    if True:
        alignSR2LC(folderName, mummerLink)
    breakLC(folderName)
        
'''
        
Subroutine 2: Obtain SC_n
    Input : SC.fasta, LC_n.fasta, LR.fasta
    Output : SC_n.fasta
    Algorithm : 
    
    0) SC_n = SC
    1) Find redundant contigs @ SC_n to drop
    2) a) Align LR to SC_n
       b) Break SC_n
    3) Do (1) again
    4) Return SC_n 
    
    

    Remark : criterion to break == long read splitting

Command to run :  python -m  srcRefactor.misassemblyFixerLib.merger
from ..repeatPhaserLib.finisherSCCoreLib import IORobot

'''

def findRedundantSC(folderName, mummerLink):
    nonRedundantResolver.removeRedundantRefvsQuery(folderName, mummerLink,  "LC_n.fasta" , "SC_n.fasta", "SC_n_tmp")

def removeNeighbor(oldList, gapBtwBk ):
    oldList.sort()
    modifiedList = [] 
    oldItm = -2*gapBtwBk
    
    for eachitem in oldList:
        if eachitem[1] > oldItm +gapBtwBk:
            modifiedList.append(eachitem)
            oldItm = eachitem[1] 
        
    return modifiedList 


    
def findPartner(seedList, oldList, gapLen):
    partnerList = []
    
    seedList.sort()
    
    lList = []
    rList = []
    
    for eachitem in oldList:
        if eachitem[2] =='L':
            lList.append(eachitem)
        elif eachitem[2] == 'R':
            rList.append(eachitem)
    
    lList.sort()
    rList.sort()
    
    for eachitem in seedList:
        if eachitem[2] == 'L':
            searchitem = [eachitem[0], eachitem[1]]
            i = bisect.bisect(rList, searchitem)
            assert(0 <= rList[i][1] - eachitem[1])
            if  rList[i][1] - eachitem[1]  < gapLen:
                insertItem = [eachitem[0], eachitem[1],rList[i][1] ,eachitem[-1],rList[i][-1] ]
                partnerList.append(insertItem)
            else:
                insertItem = [eachitem[0], eachitem[1],None ,eachitem[-1],None ]
                partnerList.append(insertItem)
                
        elif eachitem[2] == 'R':
            searchitem = [eachitem[0], eachitem[1]]
            i = bisect.bisect(lList, searchitem)
            assert(0<= eachitem[1] - lList[i-1][1])
            
            if  eachitem[1] - lList[i-1][1]  < gapLen:
                insertItem = [eachitem[0], lList[i-1][1], eachitem[1], lList[i-1][-1], eachitem[-1]]
                partnerList.append(insertItem)
            else:
                insertItem = [eachitem[0], None, eachitem[1], None, eachitem[-1]]
                partnerList.append(insertItem)
        
    return partnerList

def alignLR2SC(folderName, mummerLink):
    
    dataList = alignerRobot.largeRvsQAlign(folderName, 20, mummerLink, "SC_n", "LR", "LR2SC")
    lenDic = IORobot.obtainLength(folderName, "SC_n.fasta")
    readLenDic = IORobot.obtainLength(folderName , "LR.fasta")
    isBridgedDic = {}
    
    for eachitem in lenDic:
        isBridgedDic[eachitem] = [False for i in range(lenDic[eachitem])]
    
    thres = 30 
    boundary = 5 
    toMatchThres = 1000 
    gapBtwBk = 10000 
    gapSingleRd = 2000
    
    dataList.sort(key = itemgetter(-1), reverse= True)
    
    # Fix the reverse strand
    alignerRobot.transformCoor(dataList)
    
    bkPtList = []
    
    for key, items in groupby(dataList, itemgetter(-1)):
        ck = False
        
        readName = key
        readLen= readLenDic[readName] 
        
        tmpList = []
        for eachmaploc in items:
            matchLen = eachmaploc[5] 
            tmpList.append(eachmaploc)
            
            if abs(readLen - matchLen) < thres :
                contigName = eachmaploc[-2]
                startContig, endContig = eachmaploc[0], eachmaploc[1]
                 
                for eachbase in range( startContig + boundary , endContig-boundary):
                    isBridgedDic[contigName][eachbase] = True
                ck = True 
                break
            
        if not ck :
            
            maxLlen = 0
            maxRlen = 0 
            Rseek = []
            Lseek = []
             
            for eachmaploc2 in tmpList:
                
                startrd, endrd, matchrdlen  = eachmaploc2[2], eachmaploc2[3], eachmaploc2[5]
                
                if startrd < thres/2 and matchrdlen > maxLlen:
                    maxLlen, Rseek = matchrdlen, eachmaploc2
                if endrd > readLenDic[readName] - thres/2 and matchrdlen > maxRlen:
                    maxRlen, Lseek = matchrdlen, eachmaploc2
            
            if maxLlen > toMatchThres and maxRlen > toMatchThres:
                if abs(Rseek[1] - Lseek[0]) > gapSingleRd:
                    bkPtList.append([Rseek[-2], Rseek[1], 'R', Rseek[-1] ])
                    bkPtList.append([Lseek[-2], Lseek[0], 'L', Lseek[-1] ])
     
    breakDic = {} 
    bkPtList.sort()
    
    
    # New version : bkPtList format [[ contigName, contigBkPt, L/R dir to extend, readName]]
    #               and include top and bottom
    for key, items in groupby(bkPtList, itemgetter(0)):
        for eachitem in items:
            contigName , locToBreak = eachitem[0], eachitem[1]
            if isBridgedDic[contigName][locToBreak] == True:
                if contigName in breakDic:
                    breakDic[contigName] += [eachitem]
                else:
                    breakDic[contigName] = [eachitem] 
    
    ### remove close neighbor and  add top and bottom
    
    newBreakDic= {}
    for eachitem in breakDic:
        seedList = removeNeighbor(breakDic[eachitem], gapBtwBk )
        partnerList = findPartner(seedList, breakDic[eachitem],gapBtwBk )
        
        startItem, endItem = [eachitem, 0, 0, None, None], [eachitem, lenDic[eachitem] ,lenDic[eachitem], None, None]
        
        newBreakDic[eachitem] = [startItem]
        newBreakDic[eachitem] += partnerList
        newBreakDic[eachitem] +=  [endItem]

    # New version : breakDic format : 
    # {'seg1' : [ [contigName, bkPtL, bkPtR, readNameL, readNameR] ... ], [[],[], ... ]}
    
    with open(folderName + "SC_n_breakList.json", 'w') as outfile:
        json.dump(newBreakDic, outfile)
    
def breakSC(folderName, mummerLink):
    # New version : breakDic format : 
    # {'seg1' : [ [contigName, bkPtL, bkPtR, readNameL, readNameR] ... ], [[],[], ... ]}
    
    json_data = open(folderName + "SC_n_breakList.json", 'r')
    breakDic = json.load(json_data)
    
    contigList = IORobot.loadContigsFromFile(folderName, "SC_n.fasta")
    
    brokenContigList = []
    
    for eachitem in contigList:
        if eachitem in breakDic:
            #contigBreakDown = breakAcBkPts(contigList[eachitem], breakDic[eachitem])
            
            contigBreakDown = breakAcBkPtsTwoSided(contigList[eachitem], breakDic[eachitem],folderName, mummerLink)
            
            brokenContigList = brokenContigList + contigBreakDown
        else:
            brokenContigList = contigList[eachitem]
            
    IORobot.writeSegOut(brokenContigList, folderName, "SC_n.fasta")
     

def fixSCMisassembly(folderName , mummerLink):
    print "fixSCMisassembly"
    
    command = "cp "+ folderName + "SC.fasta "+ folderName + "SC_n.fasta"
    os.system(command)

    if True:
        findRedundantSC(folderName, mummerLink)
    
        alignLR2SC(folderName, mummerLink)
    
    breakSC(folderName, mummerLink)
    
    findRedundantSC(folderName, mummerLink)
    
        


