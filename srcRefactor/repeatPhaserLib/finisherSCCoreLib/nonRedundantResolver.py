import os
import alignerRobot
import IORobot
import houseKeeper

# ## 0) Preprocess by removing embedded contigs (I: contigs.fasta ; O : noEmbed.fasta)

def removeEmbedded(folderName , mummerLink):
    print "removeEmbedded"
    thres = 10
    os.system("sed -e 's/|//g' " + folderName + "contigs.fasta  > " + folderName + "contigs2.fasta")

    os.system("cp " + folderName + "contigs2.fasta " + folderName + "contigs.fasta") 

    if True:
        alignerRobot.useMummerAlignBatch(mummerLink, folderName, [["self", "contigs.fasta", "contigs.fasta", ""]], houseKeeper.globalParallel )
        # alignerRobot.useMummerAlign(mummerLink, folderName, "self", "contigs.fasta", "contigs.fasta")
        # outputName, referenceName, queryName, specialName
    
    dataList = alignerRobot.extractMumData(folderName, "selfOut")
    
    dataList = alignerRobot.transformCoor(dataList)
    
    lenDic = IORobot.obtainLength(folderName, 'contigs.fasta')
    
    removeList = []
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
                
    nameList = obtainComplement(lenDic, removeList)
    
    IORobot.putListToFileO(folderName, "contigs.fasta", "noEmbed", nameList)


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
    
    

    
    
