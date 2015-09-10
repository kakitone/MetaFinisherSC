from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver
from ..repeatPhaserLib.finisherSCCoreLib import alignerRobot
from ..repeatPhaserLib.finisherSCCoreLib import IORobot

from itertools import groupby
from operator import itemgetter
import merger
import os 
import json


def fixAdaptorSkip(folderName, mummerLink, inputFileName, outputFileName):
    print "fixAdaptorSkip"
    thres = 5000
    delta = 5000 
    largeThres = 5000
    
    if True:
        outputName, referenceName, queryName =  "adaptorCk", inputFileName, inputFileName
        command = mummerLink + "nucmer  -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
        os.system(command)
    
        alignerRobot.showCoorMummer(False, mummerLink, folderName, outputName, "")

    dataList = alignerRobot.extractMumData(folderName, "adaptorCkOut")
    lenDic = IORobot.obtainLength(folderName, inputFileName)
    
    dataL = len(dataList[0])
    dataList.sort(key = itemgetter(dataL-2,dataL-1))
    
    breakDic = {}
    for key, items in groupby(dataList, itemgetter(dataL-2,dataL-1)):
        #print key, lenDic[key[0]]
        if key[0] == key[1] :
            print key, lenDic[key[0]]
            contigName = key[0]
            lcontig = lenDic[contigName]
            for eachitem in items:
                print eachitem
                start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
                if start2 > end2 and (end1 - start1) > largeThres:
                    if start1 < thres:
                        if 0 <= end2 - end1 < delta: 
                            breakDic[contigName] = [0, end1, end2, lcontig] 
                            break
                        elif end1 >  lcontig - delta :
                            breakDic[contigName] = [0, end1/2-1, end1/2 +1, lcontig]
                            break 
                    elif end1 > lenDic[contigName]-thres:
                        if 0 <=start1 - start2 < delta:
                            breakDic[contigName] = [0, start2,  start1, lcontig] 
                            break
                        elif start1 < delta:
                            doublelen = end1 - start1
                            breakDic[contigName] = [0,  lcontig-doublelen/2-1, lcontig-doublelen/2 +1 ,lcontig]
                            break
                    elif start2 > lenDic[contigName] - thres:
                        start1, end1, start2, end2 = eachitem[3], eachitem[2], eachitem[1], eachitem[0]
                        #print "here"
                        print eachitem
                        print start1, end1, start2, end2
                        if 0 <=start1 - start2 < delta:
                            print "here1"
                            breakDic[contigName] = [0, start2,  start1, lcontig]
                            break
                        elif start1 < delta:
                            doublelen = end1 - start1
                            print "here2"
                            breakDic[contigName] = [0,  lcontig-doublelen/2-1, lcontig-doublelen/2 +1 ,lcontig]
                            break

                    elif end2 < thres: 
                        start1, end1, start2, end2 = eachitem[3], eachitem[2], eachitem[1], eachitem[0]    
                        if 0 <= end2 - end1 < delta:
                            print "kkhere"
                            breakDic[contigName] = [0, end1, end2, lcontig]
                            break
                        elif end1 >  lcontig - delta :
                            print "kkhere"
                            breakDic[contigName] = [0, end1/2-1, end1/2 +1, lcontig]
                            break


    contigDic = IORobot.loadContigsFromFile(folderName, inputFileName)
    returnContigList = []
    print breakDic, len(breakDic)
    for eachitem in lenDic:
        if eachitem in breakDic:
            returnList = merger.breakAcBkPts(contigDic[eachitem], breakDic[eachitem])
            returnList.pop(1) 
            returnContigList = returnContigList + returnList
        else:
            returnContigList = returnContigList + [contigDic[eachitem]]
                            
    IORobot.writeSegOut(returnContigList, folderName, "clearedAdaptor.fasta")
    
    ### log data
    adaptorSkippedLogDic = {}
    for eachitem in breakDic:
        adaptorSkippedLogDic[eachitem] = [0, breakDic[eachitem][1], breakDic[eachitem][-1]]
        
    with open(folderName + 'adaptorSkippedLogDic.json', 'w') as f:
        json.dump(adaptorSkippedLogDic, f)
    ###
    '''
    os.system("cp " + folderName + "contigs.fasta " + folderName + "contigsbackup.fasta")
    os.system("cp " + folderName + "clearedAdaptor_"+inputFileName + " " + \
              folderName + "contigs.fasta")
    
    '''

    if True:
        nonRedundantResolver.removeRedundantWithFile(folderName , mummerLink, "clearedAdaptor", "clearedAdaptorRedundant",outputFileName )


