from itertools import groupby
import os 

abunGlobalAvoidrefine = False
abunGlobalReadSearchDepth = 1
abunGlobalRRDisable = False

def replaceFiles( folderName, replacedName) :
    commandList = []
    commandList.append("cp " + folderName + "improved3.fasta " + folderName + "improved3_backup.fasta")
    commandList.append("cp " + folderName + replacedName + " "+folderName + "improved3.fasta")
    
    for eachcommand in commandList:
        print eachcommand
        os.system(eachcommand)

def parseEdgeNameToID(name, mytype):

    if mytype == 'C':
        dataInfo = name[6:].split('_')
    elif mytype == 'R':
        dataInfo = name[4:].split('_')
        
    contigNum = int(dataInfo[0])
    
    if dataInfo[1] == 'p':
        id = contigNum * 2
    elif dataInfo[1] == 'd':
        id = contigNum * 2 + 1
        
    return id 


def getDistinct(myList):
    newList = [] 
    myList.sort()
    
    for key, items in groupby(myList):
        newList.append(key)
    
    return newList


def filterData(dataList, lenDic):
    newDataList = []
    for eachitem in dataList:
        if headTailMatch(eachitem, lenDic):
            newDataList.append(eachitem)

    return newDataList


def headTailMatch(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    l1, l2 = lenDic[eachitem[-2]], lenDic[eachitem[-1]]
    thres = 40 
    
    
    diffContig , forwardStrand, headtailoverlap = False, False , False
    
    if eachitem[-2] == eachitem[-1]:
        diffContig = False
    else:
        diffContig = True

    if start2 > end2 : 
        forwardStrand = False
    else:
        forwardStrand = True
    
    
    if (start1 <= thres and end2 >= l2 - thres) or (end1 >= l1 - thres and start2 <= thres): 
        headtailoverlap = True
    else:
        headtailoverlap = False
        
    
    if diffContig and forwardStrand and headtailoverlap:
        return True
    else:
        return False
    
    
    return True

