import os
import houseKeeper
import graphLib
import alignerRobot

# ## 0) Preprocess by removing embedded contigs (I: contigs.fasta ; O : noEmbed.fasta)

def obtainLength(folderName, fileName):
    f = open(folderName + fileName, 'r')
    tmp = f.readline().rstrip()
    lenDic = {}
    tmplen = 0 
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if tmplen != 0:
                lenDic[tmpName] = tmplen
                tmplen = 0
            tmpName = tmp[1:]
        else:
            tmplen += len(tmp)
        tmp = f.readline().rstrip()
        
    lenDic[tmpName] = tmplen
    
    f.close()
    
    return lenDic
    
def findContigLength(folderName, option):
    
    if option == "contigs":
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "smaller_contigs_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
        
    else:
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "improved_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
    
    
        f = open(folderName + "relatedReads_Double.fasta", 'r')
        
    
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) > 0:
    
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()

    return contigLength

def putListToFileO(folderName, sourceFileName, targetFileName, myList):
    f = open(folderName + targetFileName + ".txt", 'w')
    for eachitem in myList:
        f.write(eachitem + '\n')
        
    f.close()
    
    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + folderName + targetFileName + ".txt " + folderName + sourceFileName + " > " + folderName + targetFileName + ".fasta"
    os.system(command)

def writeToFile_Double1(folderName, fileName1, fileName2, option="contig"):

    f2 = open(folderName + fileName2, 'w')
    fOriginal = open(folderName + fileName1, 'r')
    
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) > 0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead + tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)
    
    print "len(readSet)", len(readSet)
    
    fOriginal.close()
    
    if option == "contig":
        header = ">Contig"
    else:
        header = ">Read"
    for eachcontig, dum in zip(readSet, range(len(readSet))):
        f2.write(header + str(dum) + "_p\n")
        f2.write(eachcontig + '\n') 
        f2.write(header + str(dum) + "_d\n")
        f2.write(houseKeeper.reverseComplement(eachcontig) + '\n')
        
    f2.close()

def loadContigsFromFile(folderName, fileName):
    f = open(folderName + fileName, 'r')
    tmp = f.readline().rstrip()
    dataDic = {}
    tmpSeq = ""
    tmpName = ""
    
    while len(tmp) > 0:
        
        if tmp[0] == '>':
            if len(tmpSeq) != 0:
                dataDic[tmpName] = tmpSeq
                tmpSeq = ""
            tmpName = tmp[1:]
        else:
            tmpSeq += tmp
        tmp = f.readline().rstrip()
        
    dataDic[tmpName] = tmpSeq
    

    f.close()
    return dataDic

def writeToFile(f2, runningIndex, seq):
    f2.write(">Seg_" + str(runningIndex))
    f2.write('\n')
    f2.write(seq)
    f2.write('\n')

# ## 5) Read the contigs out (I: startList, graphNodes, ; O:improved.fasta, openZo`ne.txt)

def useAlignToGetLen(eachnode, i , nameDic, orientation, myContigsDic, readNum, folderName , mummerLink):
    indexToAddNext = eachnode.nodeIndexList[i+1]
    readNumNext = indexToAddNext / 2
    orientationNext = indexToAddNext % 2 
    
    if len(nameDic) > 0:
        orientationNext = nameDic[indexToAddNext]%2
        readNumNext = nameDic[indexToAddNext]/2
    
    if orientation == 0:
        leftSeg= myContigsDic['Contig' + str(readNum) + '_' + 'p']
    else:
        leftSeg = myContigsDic['Contig' + str(readNum) + '_' + 'd']
    
    if orientationNext == 0 :
        rightSeg = myContigsDic['Contig' + str(readNumNext) + '_' + 'p']
    else:
        rightSeg = myContigsDic['Contig' + str(readNumNext) + '_' + 'd']
        
    
    overlapInfo= align(leftSeg, rightSeg, folderName, mummerLink)
    overlapLen = overlapInfo[0]
    
    return overlapLen

def readContigOut(folderName, mummerLink, graphFileName, contigFile, outContigFile, outOpenList, nameDic={}):
    
    print "readContigOut"
     
    G = graphLib.seqGraph(0)
    G.loadFromFile(folderName, graphFileName)
    G.findStartEndList()

    myContigsDic = loadContigsFromFile(folderName, contigFile)
    
    contigUsed = [False for i in range(len(G.graphNodesList) / 2)]
     
    seqToPrint = []
    openList = []
    
    noForRevMismatch = True
    
    print len(G.graphNodesList)
    for eachnode in G.graphNodesList:
        print eachnode.nodeIndexList
        if len(eachnode.nodeIndexList) > 0:
            tmpSeq = ""
            # ## debug consistency of t/f
            ckList = []
            for dummy in eachnode.nodeIndexList:
                indexToAdd = dummy
                readNum = indexToAdd / 2
                ckList.append(contigUsed[readNum])
                
            if (len(ckList) > 0 and not all(ckList) and any(ckList)):
                noForRevMismatch = False
            
            # ## end debug 
            if contigUsed[eachnode.nodeIndexList[0]/2] == False:
                contigUsed[eachnode.nodeIndexList[0]/2] = True
                contigUsed[eachnode.nodeIndexList[-1]/2] = True
                
                for i in range(len(eachnode.nodeIndexList)):
                    
                    indexToAdd = eachnode.nodeIndexList[i]
                    readNum = indexToAdd / 2
                    orientation = indexToAdd % 2 
                        
                    #print nameDic[indexToAdd]
                    if len(nameDic) > 0:
                        orientation = nameDic[indexToAdd]%2
                        readNum = nameDic[indexToAdd]/2
                    
                    #print readNum
                    if i != len(eachnode.nodeIndexList) - 1:
                          
                        overlapLenOld = eachnode.overlapList[i]
                        
                        # Can we hijack here for the overlap Length... seems like minimal changes   
                        overlapLen =  useAlignToGetLen(eachnode, i , nameDic, orientation, myContigsDic, readNum, folderName , mummerLink)
                        # End Hijacking 
                        print overlapLen, overlapLenOld
                        
                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'p'][0:-overlapLen]
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'd'][0:-overlapLen]
                    else:
                        if orientation == 0:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'p']
                        else:
                            tmpSeq = tmpSeq + myContigsDic['Contig' + str(readNum) + '_' + 'd']
                        
                    
            if len(tmpSeq) > 0:
                if eachnode.nodeIndex in G.myStartList:
                    openList.append('Segkk' + str(len(seqToPrint)) + ',noprev')
                if eachnode.nodeIndex in G.myEndList:
                    openList.append('Segkk' + str(len(seqToPrint)) + ',nonext')
                

                seqToPrint.append(tmpSeq)


    print "No forward/reverse mismatch ?", noForRevMismatch
    fImproved = open(folderName + outContigFile, 'w')
    for eachcontig, dummyIndex in zip(seqToPrint, range(len(seqToPrint))):
        print len(eachcontig)
        fImproved.write(">Segkk" + str(dummyIndex) + '\n')
        fImproved.write(eachcontig + '\n')
         
    fImproved.close()
    
    print "All contigs used? ", all(contigUsed)
    print "NContig", len(seqToPrint)

    f = open(folderName + outOpenList, 'w')
    f.write(str(len(seqToPrint)) + '\n')
    for eachitem in openList:
        f.write(str(eachitem) + str('\n'))
    f.close()


def obtainLinkInfo(folderName, mummerLink, inputFile, mummerFile):
    # minLen = 200
    minLen = 400
    
    if houseKeeper.globalRelaxThres == False:
        thres = 5
    elif houseKeeper.globalRelaxThres == True:
        thres = 10
    
    
    writeToFile_Double1(folderName, inputFile + ".fasta", inputFile + "_Double.fasta", "contig")
    
    fmyFile = open(folderName + inputFile + "_Double.fasta", 'r')
    fSmaller = open(folderName + inputFile + "_contigs_Double.fasta", 'w')

    tmp = fmyFile.readline().rstrip()
    maxSize = 50000

    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp + '\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize)] 
            countComp = len(component)
            fSmaller.write(component)
            
            component = tmp[max(0, len(tmp) - maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            

            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()
    
    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, [[mummerFile, inputFile + "_contigs_Double.fasta", inputFile + "_contigs_Double.fasta", ""]], houseKeeper.globalParallel )
        
        # alignerRobot.useMummerAlign(mummerLink, folderName, mummerFile, inputFile + "_contigs_Double.fasta", inputFile + "_contigs_Double.fasta")
        
        
    lengthDic = obtainLength(folderName, inputFile + "_contigs_Double.fasta") 
    
    dataSetRaw = alignerRobot.extractMumData(folderName, mummerFile + "Out")
    
    # ## Format [ helperStart, helperEnd , readStart, readEnd,matchLen1,matchLen2,percentMatch,helperName,readName]
    
    
    dataSet = []
    
    for eachitem in dataSetRaw: 
        helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName = eachitem 
        
        detailHelper = helperName.split('_')
        detailRead = readName.split('_')
        

        if detailHelper[0] != detailRead[0] and  helperName != readName and max(matchLen1, matchLen2) > minLen and readStart < readEnd  and min(helperStart, readStart) < thres and min(lengthDic[helperName] - helperEnd, lengthDic[readName] - readEnd) + 1 < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if helperStart < thres:
                
                dataSet.append((max(matchLen1, matchLen2), readName, helperName))
    
    dataSet.sort(reverse=True)
    
    numberOfContig = len(lengthDic)
    
    return numberOfContig, dataSet

def truncateEndOfContigs(folderName, filenameIn, filenameOut, maxSize, lengthDic):    
    
    '''
    fmyFile = open(folderName + "improved_Double.fasta", 'r')
    fSmaller = open(folderName + "smaller_improvedContig.fasta", 'w')
    maxSize = 25000
    truncateEndOfContigs(folderName, "improved_Double.fasta", "smaller_improvedContig.fasta", 25000, lengthDic)
    
    '''
    fmyFile = open(folderName + filenameIn, 'r')
    fSmaller = open(folderName + filenameOut, 'w')
    
    tmp = fmyFile.readline().rstrip()
    
    myName = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            fSmaller.write(tmp + '\n')
            myName = tmp[1:]
        else:
            component = tmp[0:min(len(tmp), maxSize)] 
            countComp = len(component)
            fSmaller.write(component)
            
            component = tmp[max(0, len(tmp) - maxSize):len(tmp)]
            fSmaller.write(component)
            countComp = countComp + len(component)
            
            lengthDic[myName] = countComp 
            print "DebugName", myName, countComp
            fSmaller.write('\n')

        tmp = fmyFile.readline().rstrip()

    fSmaller.close()
    fmyFile.close()

def obtainLinkInfoReadContig(dummyI, mummerLink, folderName,thres, lengthDic, K):
    dataSet = []
    indexOfMum = ""
    if dummyI < 10:
        indexOfMum = "0" + str(dummyI)
    else:
        indexOfMum = str(dummyI)

        '''
        command = mummerLink + "nucmer --maxmatch --simplify -p " + folderName + "outRefine " + folderName + "smaller_improvedContig.fasta " + "relatedReads_Double.part-" + indexOfMum + ".fasta"
        os.system(command)
        
        command = mummerLink + "show-coords -r " + folderName + "outRefine.delta > " + folderName + "fromMumRefine" + indexOfMum
        os.system(command)
        '''
        
        
    f = open(folderName + "fromMumRefine" + indexOfMum, 'r')
    
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr = info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        matchLenArr = info[2].split()
       
    
        matchLen = int(matchLenArr[1])    
        contigStart, contigEnd = int(firstArr[0]), int(firstArr[1])
        readStart, readEnd = int(filterArr[0]) , int(filterArr[1])
        
        contigName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
                
        if readStart < readEnd and matchLen > K and min(contigStart, readStart) < thres and min(lengthDic[contigName] - contigEnd , lengthDic[readName] - readEnd) + 1 < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if contigStart < thres:
                dataSet.append((readName, contigName, 'L', matchLen))
            
            if lengthDic[contigName] - contigEnd + 1 < thres :
                dataSet.append((readName, contigName, 'R', matchLen))
        
        tmp = f.readline()
        
    f.close()  
    
    return dataSet

### read from overlap

def writeSegOut(oldCtgList, folderName, fileout):
    ctgList = []
    thres = 10 
    
    for eachitem in oldCtgList:
        if len(eachitem) >= thres:
            ctgList.append(eachitem)
            
    f = open(folderName + fileout, 'w')
    
    for i in range(len(ctgList)):
        f.write(">Segkk" + str(i) +'\n')
        f.write(ctgList[i])
        f.write("\n")
    
    f.close()
    
def checkIncluded(tmp, markedList):
    isIncluded  = False
    for i in tmp:
        if markedList[i/2] == True:
            isIncluded = True
    return isIncluded

def align(leftSeg, rightSeg, folderName, mummerLink):
    return alignWithName(leftSeg, rightSeg, folderName, mummerLink, "overlap")    

def alignWithName(leftSeg, rightSeg, folderName, mummerLink, nameOfOut):
    overlap = [0, 0 ] 
    lLen = 0
    f = open(folderName + nameOfOut+"leftSeg.fasta", 'w')
    f.write(">SegL\n")
    
    if len(leftSeg) < 50000:
        f.write(leftSeg)
        lLen = len(leftSeg)
    else:
        f.write(leftSeg[-50000:])
        lLen = 50000
    f.close()
    
    rLen = 0
    f = open(folderName + nameOfOut+"rightSeg.fasta", 'w')
    f.write(">SegR\n")
    if len(rightSeg) < 50000:
        f.write(rightSeg)
        rLen  = len(rightSeg)
    else:
        f.write(rightSeg[0:50000])
        rLen = 50000
        
    f.close()
    
    
    #alignerRobot.useMummerAlign(mummerLink, folderName, "overlap", "leftSeg.fasta", "rightSeg.fasta", False)
    alignerRobot.useMummerAlign(mummerLink, folderName, nameOfOut, nameOfOut+"leftSeg.fasta", nameOfOut+"rightSeg.fasta", specialForRaw = False, specialName = "", refinedVersion= True)
    dataList =  alignerRobot.extractMumData(folderName , nameOfOut+"Out")
    
    thres = 10
    
    
    if len(dataList) == 0:
        overlap = [0, 0 ]
    else:
        myMax = [0, 0]
        
        for eachitem in dataList:
            if eachitem[1] > lLen - thres and eachitem[2] < thres:
                if eachitem[5] > myMax[1]:
                    myMax[0] = eachitem[4]
                    myMax[1] = eachitem[5]
        
        overlap = myMax 
    
    return overlap 
       
def joinSeg(tmp, folderName, segLookUp, mummerLink, gapContentLookUpDic):
    
    tmpList = []
    
    ctg = ""
    
    for i in range(len(tmp)-1):
        
        begin, end = tmp[i], tmp[i+1]
        leftEnd, rightStart, middleContent = gapContentLookUpDic[str(begin) + "_" + str(end)]
        ctg = ctg + segLookUp[begin][0:leftEnd] + middleContent
        
    ctg = ctg + segLookUp[tmp[-1]]
        
    tmpList.append(ctg)
        
    return tmpList

def readContigsFromFile(folderName, filename):
    segLookUp = [] 
    f = open(folderName + filename, 'r')
    tmp = f.readline().rstrip()
    tmpStr = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpStr) > 0:
                segLookUp.append(tmpStr)
                tmpStr = ""
        else:
            tmpStr = tmpStr + tmp
        tmp = f.readline().rstrip()
     
    if len(tmpStr) > 0: 
        segLookUp.append(tmpStr)
        tmpStr = ""
    
    f.close()
    return segLookUp

def fastaContigNameIndexConversion(folderName, filename):
    i =0 
    name2Index = {}
    index2Name ={}
    
    f = open(folderName + filename, 'r')
    tmp = f.readline().rstrip()
    
    while len(tmp) > 0:
        if tmp[0] == '>':
            name = tmp[1:].rstrip()
            name2Index[name] = i
            index2Name[i] = name 
            i = i+1
        tmp = f.readline().rstrip()
    
    f.close()
    
    return name2Index, index2Name
    
def extractGraphToContigs(G, folderName, mummerLink, fileout, filein, gapContentLookUpDic, mapDummyToRealDic= {}):
    
    N1 = len(G.graphNodesList) - len(mapDummyToRealDic)
    markedList = [False for i in range(N1/2)]
    
    segLookUp = readContigsFromFile(folderName, filein)
    
    ctgList = []
    for eachnode in G.graphNodesList:
        if eachnode.nodeIndex < N1:
            if len(eachnode.nodeIndexList) > 0: 
                myNodeIndexList = eachnode.nodeIndexList 
                print myNodeIndexList
                toCkIncludedList = []
                tmp = []
                
                for kk in myNodeIndexList:
                    if kk < N1:
                        toCkIncludedList.append(kk)
                
                for kk in myNodeIndexList:
                    if kk < N1:
                        tmp.append(kk)
                    else:
                        # tmp.append(mapDummyToRealDic[str(kk-N1)])
                        tmp = tmp + mapDummyToRealDic[str(kk-N1)][1]
                
                
                isIncluded  = checkIncluded(toCkIncludedList, markedList)
                
                #print isIncluded, toCkIncludedList, markedList
                if not isIncluded :
                    for eachitem in toCkIncludedList: 
                        markedList[eachitem/2] = True

                    ctgtmpList = joinSeg(tmp, folderName, segLookUp, mummerLink, gapContentLookUpDic)
                    ctgList = ctgList + ctgtmpList
    
    missingSegments = []
    for i in range(len(markedList)):
        if markedList[i] == False: 
            missingSegments.append(segLookUp[2*i])
            #print i
    #assert(False)
    print "len(missingSegments): " , len(missingSegments)
    
    ctgList = ctgList + missingSegments
    
    writeSegOut(ctgList, folderName, fileout)
            
def myRead(folderName, fileName,readName) :
    f = open(folderName + fileName, 'r')
    tmp = f.readline().rstrip()
    returnStr = ""
    ck = False
    while len(tmp) > 0:
        
        if tmp[0] == '>' and len(returnStr) == 0 and tmp[1:] == readName:
            ck = True
        elif tmp[0] == '>' and len(returnStr) >0 :
            ck = False
        elif tmp[0] != '>':
            if ck :
                returnStr = returnStr + tmp 
        
        tmp = f.readline().rstrip()
        
    f.close()
    
    return returnStr

def pathListToSeqListTransform(pathList, contigList, readList, mummerPath, folderName):
    pathSeqList  = []

    for p in pathList:
        pSeq = alignAndMerge(p, contigList, readList, mummerPath, folderName)
        pathSeqList.append(pSeq)

    return pathSeqList

def alignAndMerge(p, contigList, readList, mummerPath, folderName):
    pSeq = []

    N1 = len(contigList)

    pSeg = contigList[p[0]]

    for i in range(len(p)-1):
        if p[i] < N1:
            leftSeg = contigList[p[i]]
        else:
            leftSeg = readList[p[i]- N1]

        if p[i+1] < N1:
            rightSeg = contigList[p[i+1]]
        else:
            rightSeg = readList[p[i+1]- N1]

        overlap = align(leftSeg,rightSeg,folderName,mummerPath)

        rightStart = overlap[1]
        pSeg = pSeg + rightSeg[rightStart:]

    return pSeg


