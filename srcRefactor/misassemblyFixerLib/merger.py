from ..repeatPhaserLib.finisherSCCoreLib import IORobot
from ..repeatPhaserLib.finisherSCCoreLib import alignerRobot
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper 
import bisect
import adaptorFix

mergerGlobalLCReads = "SR"

class fixerRobot:
    def __init__(self):
        print "fixerRobot"
        self.toRunAdaptor = True
        self.toRunNoEmbedEnd = True
        self.toRunAggressive = False
        self.sdMult = 5
        self.tuneParaOnly = False


    def loadData(self, initial_data):
        canLoad = True
        for key in initial_data:
            if hasattr(self, key):
                if initial_data[key] =='True' :
                    setattr(self, key, True)
                elif initial_data[key] == 'False':
                    setattr(self, key, False)
                else:
                    setattr(self, key, float(initial_data[key]))
            else:
                canLoad = False
        return canLoad


mergerGlobalFixerRobot = fixerRobot()


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
    commandList = []
    command = "cat "+ folderName + "LC_n.fasta "+ folderName + "SC_n.fasta > "
    command = command + folderName + "contigs.fasta"  
    
    commandList.append(command)
    
    command = "perl -pe 's/>[^\$]*$/\">Seg\" . ++$n .\"\n\"/ge' "+folderName+"contigs.fasta > "+folderName+"newContigs.fasta "
    commandList.append(command)
    
    command = "cp " +folderName+"newContigs.fasta  "+folderName+"contigs.fasta "
    commandList.append(command)
    
    if True:  
        for eachcomm in commandList:
            
            os.system(eachcomm)
            
        
def mergeContigs(folderName, mummerLink, inputLCname):
    print "mergeContigs" 
    
    fixLCMisassembly(folderName, mummerLink, inputLCname)
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
    
    

def alignSR2LC(folderName, mummerLink, incontigName):
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
    LCLenDic = IORobot.obtainLength(folderName, incontigName+".fasta")
    coveragePerContigsDic = {}
    for eachitem in LCLenDic : 
        coveragePerContigsDic[eachitem] = [ 0 for i in range(LCLenDic[eachitem]) ]
    
    random.seed(0)
    outputName, referenceName, queryName = "SR2"+incontigName+"Align" , incontigName, "SR" 
    
    if True:
        dataList = alignerRobot.largeRvsQAlign(folderName, 20, mummerLink, referenceName, queryName, outputName)

    
    dataList.sort(key = itemgetter(-1))
    LCLenDic = IORobot.obtainLength(folderName, incontigName+ ".fasta")
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
    if mergerGlobalLCReads == "SR":
        thresValue = 13*sd
    elif mergerGlobalLCReads == "LR":
        thresValue = mergerGlobalFixerRobot.sdMult*sd
    
    n = len(tDiffVec)
    deletedList = [False for i in range(n)]
    searchList = [[abs(tDiffVec[i]), i ] for i in range(n)]
    searchList.sort(reverse = True)
    print searchList[0:10]
    
    for eachitem in searchList:
        value, index = eachitem[0] , eachitem[1]
        #print value, index
        if value < thresValue:
            break
        else:
            if deletedList[index] == True:
                continue
            else:
                outLiners.append(index)
                for j in range(max(index - T,0),  min(index+ T+1,n)):
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
   
   
    if mergerGlobalLCReads == "LR":
        return x
    elif mergerGlobalLCReads == "SR": 
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


    

def breakAcBkPtsTwoSided(contig, modifiedOutlinersOld, folderName, mummerLink):
    contigBreakDown = []
    
    
    modifiedOutliners = [] 
    for eachitem in modifiedOutlinersOld:
        tmpList = []
        if eachitem[1] != None:
            tmpList = tmpList + eachitem + [eachitem[1]]
        elif eachitem[2] != None:  
            tmpList = tmpList + eachitem + [eachitem[2]]
        modifiedOutliners.append(tmpList)
    
    modifiedOutliners.sort(key = itemgetter(-1))
    
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
                        
            if overlap[0] >0:
                leftpart = leftSeg[0:len(leftSeg)-overlap[0]]
                print "overlap:", overlap
            else:
                leftSeg = houseKeeper.reverseComplement(leftSeg)
                overlap = IORobot.align(leftSeg, rightSeg, folderName, mummerLink)
                
                if overlap[0] > 0:
                    leftpart = leftSeg[0:len(leftSeg)-overlap[0]]
                    print "overlap:", overlap
                else:
                    leftpart = ""
                    print "overlap error ? :", overlap 
                
          
        if rightRead != None:
            leftSeg, rightSeg = contig[begin:end], IORobot.myRead(folderName, "LR.fasta",rightRead ) 
            overlap = IORobot.align(leftSeg, rightSeg, folderName, mummerLink)

            if overlap[1] > 0:
                rightpart = rightSeg[overlap[1]:]
                print "overlap:", overlap
            else:
                rightSeg =  houseKeeper.reverseComplement(rightSeg)
                overlap = IORobot.align(leftSeg, rightSeg , folderName, mummerLink)
                
                if overlap[1] > 0:
                    rightpart = rightSeg[overlap[1]:]
                    print "overlap:", overlap
                else:
                    rightpart = ""
                    print "overlap error ? :", overlap 
            
        
        #assert(1==2)
        middlepart = contig[begin:end] 
        
        combinedSeg = leftpart + middlepart + rightpart
        print "len(combinedSeg)", len(combinedSeg)
        contigBreakDown.append(combinedSeg)
        
    return contigBreakDown 

def breakLC(folderName, inputName ):
    
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
    
    LCList = IORobot.readContigsFromFile(folderName, inputName+".fasta")
    lenDic = IORobot.obtainLength(folderName, inputName+".fasta")
    
    name2Index, index2Name = IORobot.fastaContigNameIndexConversion(folderName, inputName+".fasta")
    
    breakPtsDic = {}
    
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
            
            breakPtsDic[name] = modifiedOutliners
            
            contigList = contigList + contigBreakDown

        else:
            contigList = contigList + [LCList[contigIndex]]
            
    
    if mergerGlobalLCReads == "SR":
        IORobot.writeSegOut(contigList, folderName, "LC_n.fasta")
    
    with open(folderName + "modifiedOutliners.json", 'w') as outfile:
        json.dump(breakPtsDic, outfile)
    
    

def fixLCMisassembly(folderName, mummerLink, inputName):
    print "fixLCMisassembly"
    # Remark : if done, one can also try on LR only data set to see if there 
    # is any reduction on mis-assembly. 
    if True:
        alignSR2LC(folderName, mummerLink,inputName )
    breakLC(folderName, inputName)
        
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
    if True:
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
            # assert(0 <= rList[i][1] - eachitem[1])
            if  i< len(rList) and rList[i][1] - eachitem[1]  < gapLen:
                insertItem = [eachitem[0], eachitem[1],rList[i][1] ,eachitem[-1],rList[i][-1] ]
                partnerList.append(insertItem)
            else:
                insertItem = [eachitem[0], eachitem[1],None ,eachitem[-1],None ]
                partnerList.append(insertItem)
                
        elif eachitem[2] == 'R':
            searchitem = [eachitem[0], eachitem[1]]
            i = bisect.bisect(lList, searchitem)
            #assert(0<= eachitem[1] - lList[i-1][1])
            
            if  i> 0 and eachitem[1] - lList[i-1][1]  < gapLen:
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
    boundary = 50 
    toMatchThres = 1000 
    gapBtwBk = 10000 
    gapSingleRd = 50000
    
    dataList.sort(key = itemgetter(-1), reverse= True)
    
    # Fix the reverse strand
    alignerRobot.transformCoor(dataList)
    
    bkPtList = []
    
    print "dataList[0:10]", dataList[0:10]
    
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
                if abs(Rseek[1] - Lseek[0]) > gapSingleRd or Rseek[-1] !=Lseek[-1]:
                    bkPtList.append([Rseek[-2], Rseek[1]-1, 'R', Rseek[-1] ])
                    bkPtList.append([Lseek[-2], Lseek[0]-1, 'L', Lseek[-1] ])
     
    breakDic = {} 
    bkPtList.sort()
    
    print "bkPtList", bkPtList
    
    # New version : bkPtList format [[ contigName, contigBkPt, L/R dir to extend, readName]]
    #               and include top and bottom
    for key, items in groupby(bkPtList, itemgetter(0)):
        for eachitem in items:
            contigName , locToBreak = eachitem[0], eachitem[1]
            #if isBridgedDic[contigName][locToBreak] == False:
            if contigName in breakDic:
                breakDic[contigName] += [eachitem]
            else:
                breakDic[contigName] = [eachitem] 

    ### remove close neighbor and  add top and bottom
    
    print "breakDic", breakDic
    
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
            brokenContigList = brokenContigList + [contigList[eachitem]]
            
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
    
        
'''
Only Long reads and long contigs case:
'''

def onlyLRMiassemblyFix(folderName, mummerLink, inputName ):

    
    if not mergerGlobalFixerRobot.tuneParaOnly:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, [["self"+inputName, inputName+".fasta", inputName+".fasta", ""]], houseKeeper.globalParallel )
        
    dataList = alignerRobot.extractMumData(folderName, "self"+inputName+"Out")
    dataList = alignerRobot.transformCoor(dataList)
    lenDic = IORobot.obtainLength(folderName, inputName+'.fasta')
    matchThres = 10000
    nonMatchThres = 500
    count = 0

    newDataList= []
    for eachitem in dataList:
        name1, name2 = eachitem[-2], eachitem[-1]
        matchLen1 , matchLen2 = eachitem[4], eachitem[5]
        start1 , end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]

        # if name1!= name2  and matchLen1> matchThres:
        if name1!= name2   and  ( min(lenDic[name1] - end1, lenDic[name2] - end2 ) > nonMatchThres or min(start1, start2) > nonMatchThres ) and matchLen1> matchThres:
            print "eachitem ", eachitem, lenDic[name1], lenDic[name2]
            count = count + 1
            newDataList.append(eachitem)

    print "Count: " + str(count)

    if mergerGlobalFixerRobot.toRunAggressive == False:
        if not mergerGlobalFixerRobot.tuneParaOnly:
            alignSR2LC(folderName, mummerLink, inputName)
        breakLC(folderName, inputName)
        blkDic = getBreakPointFromDataList(folderName, newDataList, inputName)
    else:
        blkDic = breakPtGettingHack2(folderName, newDataList, inputName)


    LCList = IORobot.loadContigsFromFile(folderName, inputName+".fasta")

    contigList = []

    for eachcontig in LCList:

        if not eachcontig in blkDic:
            contigList = contigList + [ LCList[eachcontig] ]
        else:
            contigList = contigList + breakAcBkPts(LCList[eachcontig], blkDic[eachcontig])

    print "len(contigList)", len(contigList)
    IORobot.writeSegOut(contigList, folderName, "LC_n.fasta")
    
def withinBound(sep, mylist, bkpt ):
    mylist.sort()
    i = bisect.bisect(mylist, bkpt)
    ck = False

    if i>0 and abs(bkpt-mylist[i-1])< sep:
        ck = True

    if i < len(mylist) and abs(bkpt - mylist[i]) <sep :
        ck = True

    return ck


def getBreakPointFromDataList(folderName, dataList, inputLCname):
    g = 1000
    blkDic = {}
    dataList.sort(key = itemgetter(-2))
    lenDic = IORobot.obtainLength(folderName, inputLCname+ ".fasta")

    json_data = open(folderName + "modifiedOutliners.json", 'r')
    breakPtsDic = json.load(json_data)
    sep = 5000

    for key, items in groupby(dataList,itemgetter(-2) ):
        contigName = key
        newList =[]
        for eachitem in items:
            newList.append([eachitem[0], eachitem[1]])
        newList.sort()

        bktmp = [0]

        if newList[0][0] > g :
            if withinBound(sep, breakPtsDic[contigName],newList[0][0]):
                bktmp.append(newList[0][0])

        for i in range(len(newList)-1):
            if newList[i+1][0] > newList[i][1] + g:
                if withinBound(sep, breakPtsDic[contigName],newList[i+1][0]):
                    bktmp.append(newList[i+1][0])

        bktmp.append(lenDic[contigName])

        blkDic[contigName] = bktmp
        print "contigName: "+ contigName
        print "bktmp:", bktmp
        print "breakPtsDic[contigName]",breakPtsDic[contigName]

    return blkDic

def breakPtGettingHack(folderName, dataList, inputLCname):
    
    blkDic = {}
    dataList.sort(key = itemgetter(-2))
    lenDic = IORobot.obtainLength(folderName, inputLCname+".fasta") 
    g = 1000

    for key, items in groupby(dataList,itemgetter(-2) ):
        contigName = key
        newList =[]

        for eachitem in items:
            newList.append([eachitem[0], eachitem[1]])

        newList.sort()

        bktmp = [0]

        if newList[0][0] > g :
            bktmp.append(newList[0][0])

        for i in range(len(newList)-1):
            if newList[i+1][0] > newList[i][1] + g:
                bktmp.append(newList[i+1][0])

        bktmp.append(lenDic[contigName])

        blkDic[contigName] = bktmp
        print "contigName: "+ contigName
        print "bktmp:", bktmp
        #print "breakPtsDic[contigName]",breakPtsDic[contigName]

    return blkDic


def breakPtGettingHack2(folderName, dataList, inputLCname):
    
    blkDic = {}
    dataList.sort(key = itemgetter(-2))
    lenDic = IORobot.obtainLength(folderName, inputLCname+".fasta") 
    g = 1000

    for key, items in groupby(dataList,itemgetter(-2) ):
        contigName = key
        newList =[]

        for eachitem in items:
            newList.append([eachitem[0], eachitem[1]])
            
        newList.sort()
        
        bktmp = [0]

        if newList[0][0] > g :
            bktmp.append(newList[0][0])
            bktmp.append(newList[0][1])

        minThres = newList[0][1] + g

        for i in range(len(newList)-1):
            if newList[i+1][0] > minThres:
                bktmp.append(newList[i+1][0])
                bktmp.append(newList[i+1][1])

            if newList[i+1][1] + g > minThres: 
                minThres = newList[i+1][1] + g

        bktmp.append(lenDic[contigName])

        blkDic[contigName] = bktmp
        print "contigName: "+ contigName
        print "bktmp:", bktmp
        #print "breakPtsDic[contigName]",breakPtsDic[contigName]

    return blkDic


def mainFlow(newFolderName, newMummerLink):

    filterName = "LC_filtered"
    finalOutName = "mFixed"

    if not mergerGlobalFixerRobot.tuneParaOnly: 
        if mergerGlobalFixerRobot.toRunAdaptor == True:
            adaptorFix.fixAdaptorSkip(newFolderName, newMummerLink, "LC.fasta", filterName) 
        else:
            os.system("cp "+ newFolderName + "LC.fasta "+newFolderName+filterName+".fasta")


    if mergerGlobalLCReads == "SR":
        mergeContigs(newFolderName , newMummerLink, filterName)
        command = "cp "+ newFolderName + "LR.fasta "+ newFolderName + "raw_reads.fasta"
        os.system(command)
        print "Command: ",  command 
        
    elif mergerGlobalLCReads == "LR":
        
        command = "cp "+ newFolderName + "LR.fasta " + newFolderName + "SR.fasta "
        os.system(command)
        print "Command: ",  command 
        
        onlyLRMiassemblyFix(newFolderName, newMummerLink, filterName)
        
        command = "cp "+ newFolderName + "LC_n.fasta "+ newFolderName + "contigs.fasta"
        os.system(command)
        print "Command: ",  command 

        command = "cp "+ newFolderName + "LR.fasta "+ newFolderName + "raw_reads.fasta"
        os.system(command)
        print "Command: ",  command 
        

    nonRedundantResolver.removeRedundantWithFile(newFolderName , newMummerLink, "contigs", "mFixingFinalRedundantRemove", finalOutName)



