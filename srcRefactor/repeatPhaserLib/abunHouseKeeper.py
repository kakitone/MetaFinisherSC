from itertools import groupby
import os 
from finisherSCCoreLib import IORobot

abunGlobalAvoidrefine = True
abunGlobalReadSearchDepth = 0
abunGlobalRRDisable = True
abunGlobalRunPickUp = "map"

class abunSplitParameterRobot():
    def __init__(self):
        self.runGraphSurgery = True
        self.runBResolve = True
        self.runXResolve = True

        self.parameterForGraphSurgery()
        self.parameterForBResolve()
        self.parameterForXResolve()

        self.BRThres = 2
        self.AbunLower = 0.5
        self.AbunUpper = 1.95

        self.avoidrefine = abunGlobalAvoidrefine
        self.readSearchDepth = abunGlobalReadSearchDepth
        self.rrDisable = abunGlobalRRDisable

    def parameterForGraphSurgery(self):
        self.edgeThres = 1
        self.kthres = 3

        self.toRunCondenseRemove = True
        self.toRunTransitive = True
        self.toRunDoubltPtr = True


    def parameterForBResolve(self):
        self.toRunAggB = False
        self.toRunBRB = True
        self.toRunAbunB = True
        self.RThres = 5

        self.BRThresB = -1
        self.AbunLowerB = -1
        self.AbunUpperB = -1
    

    def parameterForXResolve(self):
        self.toRunAggX = False
        self.toRunBRX = True
        self.toRunAbunX = True

        self.BRThresX = -1
        self.AbunLowerX = -1
        self.AbunUpperX = -1

    def loadData(self, initial_data):
        canLoad = True
        self.avoidrefine = abunGlobalAvoidrefine
        self.readSearchDepth = abunGlobalReadSearchDepth
        self.rrDisable = abunGlobalRRDisable

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

  
abunGlobalSplitParameterRobot = abunSplitParameterRobot()


def replaceFiles( folderName, replacedName) :
    commandList = []
    commandList.append("cp " + folderName + "improved3.fasta " + folderName + "improved3_backup.fasta")
    commandList.append("cp " + folderName + "improved3_Double.fasta " + folderName + "improved3_backup.fasta")
    
    IORobot.writeToFile_Double1(folderName, replacedName[0:-6]+".fasta", replacedName[0:-6]+"_Double.fasta", "contig")
    
    commandList.append("cp " + folderName + replacedName + " "+folderName + "improved3.fasta")
    
    command = "perl -pe 's/>[^\$]*$/\">Segkk\" . $n++ .\"\n\"/ge' "+folderName+"improved3.fasta > "+folderName+"newImproved3.fasta "
    commandList.append(command)
    
    command = "cp " +folderName+"newImproved3.fasta  "+folderName+"improved3.fasta "
    commandList.append(command)


    commandList.append("cp " + folderName + replacedName[0:-6]+"_Double.fasta " + folderName + "improved3_Double.fasta")
    
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



def parseIDToName(id, mytype, numDoubleContigs):

    if mytype == 'C':
        name = "Contig"+  str(id/2)
    elif mytype == 'R':
        name =  "Read" + str((id-numDoubleContigs)/2)
    
    
    if id%2 ==0 : 
        name = name + "_p"
    else:
        name= name + "_d"

    return name 

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

def filterDataIdentical(dataList, lenDic):
    newDataList = []
    for eachitem in dataList:
        if not identicalItem(eachitem, lenDic):
            newDataList.append(eachitem)

    return newDataList

def identicalItem(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    thres = 20 
    
    endPt1 = lenDic[eachitem[-2]] - end1
    endPt2 = lenDic[eachitem[-1]] - end2
    
    if abs(start1-start2) < thres and abs(endPt1 - endPt2) < thres:
        return True
    else:
        return False

def headTailMatch(eachitem, lenDic):
    start1, end1, start2, end2 = eachitem[0], eachitem[1], eachitem[2], eachitem[3]
    l1, l2 = lenDic[eachitem[-2]], lenDic[eachitem[-1]]
    thres = 10 
    
    
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

