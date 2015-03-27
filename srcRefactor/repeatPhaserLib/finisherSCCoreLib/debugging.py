from multiprocessing import Pool
import time 
import IORobot
import alignerRobot
from itertools import groupby
from operator import itemgetter


def f(x, kk):
    print "in"
    t =0 
    for i in range(110):
        t = ((t + i) %10 * 10 )  /10 + 10
    time.sleep(1)
    
    print x + kk
    
    
def cube(x):
    return x**3



def test1():
    t= time.time()
    p = Pool(4)
    kk = 10
    results = []
    for x in range(1,7):
        results.append(p.apply_async(f, args=(x,kk)))
    output = [p.get() for p in results]
    print output
    print time.time() - t

def test2():
    pool = Pool(processes=4)
    results = [pool.apply_async(cube, args=(x,)) for x in range(1,100)]
    print results
    output = [p.get() for p in results]
    print(output)
    
    
def tmpBreakAcBkPts(contig, modifiedOutliners):
    contigBreakDown = []
    m = len(modifiedOutliners) - 2
    x = sorted(modifiedOutliners)
    
    for i in range(m+1):
        start , end = x[i], x[i+1]
        newItem = contig[start:end]
        contigBreakDown.append(newItem)
        
    return contigBreakDown


def observeOverlap(folderName):
    dataList = alignerRobot.extractMumData(folderName, "selfOut")
    dataList = alignerRobot.transformCoor(dataList)
    lenDic = IORobot.obtainLength(folderName, 'contigs.fasta')
    matchThres = 10000
    nonMatchThres = 500
    count = 0
    
    newDataList= [] 
    for eachitem in dataList:
        name1, name2 = eachitem[-2], eachitem[-1]
        matchLen1 , matchLen2 = eachitem[4], eachitem[5]
        start1 , end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
        if name1!= name2 and ( min(lenDic[name1] - end1, lenDic[name2] - end2 ) > nonMatchThres \
        or min(start1, start2) > nonMatchThres ) \
        and matchLen1> matchThres:
            print "eachitem ", eachitem, lenDic[name1], lenDic[name2]
            count = count + 1
            newDataList.append(eachitem)

    print "Count: " + str(count)
    
    blkDic = getBreakPointFromDataList(fodlerName, newDataList)
    
    LCList = IORobot.loadContigsFromFile(folderName, "contigs.fasta")
    
    contigList = []
    for eachcontig in LCList:
        if not eachcontig in blkDic:
            contigList = contigList + [ LCList[eachcontig] ]
        else:
            contigList = contigList + tmpBreakAcBkPts(LCList[eachcontig], blkDic[eachcontig])
    
    
    IORobot.writeSegOut(contigList, folderName, "breakChains.fasta")

def getBreakPointFromDataList(folderName, dataList):
    g = 1000 
    blkDic = {}
    dataList.sort(key = itemgetter(-2))
    lenDic = IORobot.obtainLength(folderName, "contigs.fasta")
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
        
    return blkDic
    
        
        
        

observeOverlap("/data/kakitone/Feb12-2015/workingMetaFinisherSC/dataFolder/")