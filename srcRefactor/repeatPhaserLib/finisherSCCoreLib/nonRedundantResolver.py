import os
import alignerRobot
import IORobot
import houseKeeper

# ## 0) Preprocess by removing embedded contigs (I: contigs.fasta ; O : noEmbed.fasta)

### Disjoint Union Data Structure
class clusterElem(object):
    def __init__(self,index):
        self.rank = 0
        self.parent = self
        self.id = index
        self.childList =[]
        
        #self.size = 1

def find(x):
    #if x != x.parent:
    #    x.parent = find(x.parent)
    #return x.parent
    if x.parent == x:
        return x
    else:
        return find(x.parent)
       
def union(x,y):
    xRoot = find(x)
    yRoot = find(y)
    
    if xRoot == yRoot:
        return 0
    
    if xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
        yRoot.childList.append(xRoot)
        
        
    elif xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
        
        
    else:
        yRoot.parent = xRoot
        xRoot.childList.append(yRoot)
        
        xRoot.rank = xRoot.rank + 1 
        
    return 1

def familyList(x):
    root = find(x)
    stack = []
    familyKmers = [root]
    
    stack.append(root)
    while (len(stack) >0 ):
        item = stack.pop(0)
        for eachsubitem in item.childList:
            familyKmers.append(eachsubitem)
            stack.append(eachsubitem)
            
    return familyKmers
                    




def removeEmbedded(folderName , mummerLink):
    print "removeEmbedded"
    removeRedundantWithFile(folderName , mummerLink, "contigs", "self", "noEmbed")


def removeRedundantWithFile(folderName , mummerLink, inputFilename, mummerTmpName, outputFileName):
    thres = 10
    os.system("sed -e 's/|//g' " + folderName + inputFilename+".fasta  > " + folderName + inputFilename+ "2.fasta")

    os.system("cp " + folderName + inputFilename+"2.fasta " + folderName + inputFilename+".fasta") 

    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, [[mummerTmpName, inputFilename+".fasta", inputFilename+".fasta", ""]], houseKeeper.globalParallel )
        # alignerRobot.useMummerAlign(mummerLink, folderName, "self", "contigs.fasta", "contigs.fasta")
        # outputName, referenceName, queryName, specialName
    
    dataList = alignerRobot.extractMumData(folderName, mummerTmpName+ "Out")
    
    dataList = alignerRobot.transformCoor(dataList)
    
    lenDic = IORobot.obtainLength(folderName, inputFilename+'.fasta')
    
    removeList = []

    shortEmbedClusterDic = {}

    for eachitem in lenDic:
        shortEmbedClusterDic[eachitem]= clusterElem(eachitem)


    for eachitem in dataList:
        match1, match2, name1, name2 = eachitem[4], eachitem[5], eachitem[7], eachitem[8]
        
        if name1 != name2:
            l1, l2 = lenDic[name1], lenDic[name2]
            
            if abs(l1 - match1) < thres and abs(l2 - match2) > thres:
                removeList.append(name1)
            elif abs(l1 - match1) > thres and abs(l2 - match2) < thres:
                removeList.append(name2)
            elif abs(l1 - match1) < thres and abs(l2 - match2) < thres:
                print "Both shortembedd", eachitem
                union(shortEmbedClusterDic[name1], shortEmbedClusterDic[name2])

        
    nameList = obtainComplement(lenDic, removeList)
    
    returnList = []

    for eachitem in nameList:
        if find(shortEmbedClusterDic[eachitem]).id == eachitem:
            returnList.append(eachitem)

    print "len(nameList), len(returnList)", len(nameList), len(returnList)

    IORobot.putListToFileO(folderName, inputFilename+".fasta", outputFileName, returnList)


def obtainComplement(lenDic, removeList):
    nameList = []
    for eachitem in lenDic:
        nameList.append(eachitem)

    print len(nameList)
    
    for eachitem in removeList:
        if eachitem in nameList:
            nameList.remove(eachitem)
    print len(nameList)
    return nameList

def removeRedundantRefvsQuery(folderName, mummerLink,  fileR , fileQ, outputFileName):
    
    thres = 10

    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, [["redundantRvsQ", fileR, fileQ, ""]], houseKeeper.globalParallel )
    
    dataList = alignerRobot.extractMumData(folderName, "redundantRvsQOut")
    lenDicR = IORobot.obtainLength(folderName, fileR)
    lenDicQ = IORobot.obtainLength(folderName, fileQ)
    
    isRedundantList= []
    
    for eachitem in dataList:
        match1, match2, name1, name2 = eachitem[4], eachitem[5], eachitem[7], eachitem[8]
        l1, l2 = lenDicR[name1] , lenDicQ[name2]
        
        if abs(l2 - match2) < thres:
            isRedundantList.append(name2)
    
    #print lenDicQ

    nonRedundantList = obtainComplement(lenDicQ, isRedundantList)
    
    print nonRedundantList
    IORobot.putListToFileO(folderName, fileQ, outputFileName, nonRedundantList)
        
    os.system("cp "+ folderName + "SC_n_tmp.fasta "+ folderName + "SC_n.fasta")
    
    

    
    
