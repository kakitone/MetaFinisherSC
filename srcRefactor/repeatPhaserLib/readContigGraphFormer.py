import os
import sys 
from itertools import groupby
from operator import itemgetter

import abunHouseKeeper

from finisherSCCoreLib import IORobot
from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import houseKeeper


def addDataToList(dataList, G, startIndex1, startIndex2, type1, type2):

    threshold = 50
    
    for eachitem in dataList:
        wt = min(eachitem[4] , eachitem[5])
        
        if eachitem[0] < threshold:

            j = abunHouseKeeper.parseEdgeNameToID(eachitem[-2], type1) + startIndex1
            i = abunHouseKeeper.parseEdgeNameToID(eachitem[-1], type2) + startIndex2
        else:
            j = abunHouseKeeper.parseEdgeNameToID(eachitem[-1], type2) + startIndex2
            i = abunHouseKeeper.parseEdgeNameToID(eachitem[-2], type1) + startIndex1
        
        
        G.insertEdge(i, j, wt) 

# ## Debug check 
def checkGraphLength(G, N1, lenDicRR):
    for eachitem in G.graphNodesList:
        for eachnext in eachitem.listOfNextNodes:
            index, wt = eachnext[0], eachnext[1]
            if index >= N1 and eachitem.nodeIndex >= N1:
                header = "Read" + str((index - N1) / 2) + "_"
                if (index - N1) % 2 == 0:
                    header = header + "p"
                else:
                    header = header + "d"
                
                header2 = "Read" + str((eachitem.nodeIndex - N1) / 2) + "_"
                    
                if (eachitem.nodeIndex - N1) % 2 == 0:
                    header2 = header2 + "p"
                else:
                    header2 = header2 + "d"
                    
                if lenDicRR[header] < wt:
                    print lenDicRR[header], lenDicRR[header2], wt
                    print header, header2
   
def alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header ):   
    #alignerRobot.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)
    numberOfFiles = 20
    bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
    command = bindir + "/finisherSCCoreLib/fasta-splitter.pl --n-parts " + str(numberOfFiles) + " " + folderName + queryFile
    os.system(command)
    
    os.system("cp *.fasta " + folderName )
    os.system("rm *.fasta ")
    
    workerList = []
    
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
       
        outputName, referenceName, queryName, specialName= header+indexOfMum, referenceFile,queryFile[0:-6]+".part-"+ indexOfMum + ".fasta" ,  header + indexOfMum
        workerList.append([outputName, referenceName, queryName, specialName])
        
    alignerRobot.useMummerAlignBatch(mummerLink, folderName, workerList, houseKeeper.globalParallel ,specialForRaw = False, refinedVersion = False)
    alignerRobot.combineMultipleCoorMum( True, mummerLink, folderName, header,header +"Out", numberOfFiles)
                    

def formReadContigStringGraph(folderName, mummerLink, contigFilename, readsetFilename, optTypeFileHeader, graphName, needAlignment=True):
    
    '''
    Input : all_associated_reads.fasta, improved3.fasta
    Output : (G) String Graph linking the reads and contigs
    Algorithm: 
        a) Form double reads and contigs                            V
        b) Mummer the data and extract dataList three times         V
        c) Use the subroutine to output a graph                     V
        d) Output the graph to a file phasing_String_graph.graph    V
    '''

    G = []

    IORobot.writeToFile_Double1(folderName, contigFilename + ".fasta", contigFilename + "_Double.fasta", "contig")
    IORobot.writeToFile_Double1(folderName, readsetFilename + ".fasta", readsetFilename + "_Double.fasta", "reads")
    
    
    header, referenceFile, queryFile = optTypeFileHeader + "CC", contigFilename + "_Double.fasta" , contigFilename + "_Double.fasta"
    if needAlignment:
        alignerRobot.useMummerAlign(mummerLink, folderName, header, referenceFile, queryFile)

    lenDicCC = IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")
    dataListCC = alignerRobot.extractMumData(folderName, header + "Out")
    dataListCC = abunHouseKeeper.filterData(dataListCC, lenDicCC)
    
    header, referenceFile, queryFile = optTypeFileHeader + "RR", readsetFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    
    
    lenDicRR = IORobot.obtainLength(folderName, readsetFilename + "_Double.fasta")
    
    if not abunHouseKeeper.abunGlobalRRDisable:
        if needAlignment:
            alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header )
    
        dataListRR = alignerRobot.extractMumData(folderName, header + "Out")
        dataListRR = abunHouseKeeper.filterData(dataListRR, lenDicRR)
        dataListRR = abunHouseKeeper.filterDataIdentical(dataListRR, lenDicRR)

    else:
        dataListRR = []
    
    header, referenceFile, queryFile = optTypeFileHeader + "CR", contigFilename + "_Double.fasta" , readsetFilename + "_Double.fasta"
    if needAlignment:
        alignerSubRoutine(folderName ,referenceFile,  queryFile, mummerLink, header )
    
    lenDicCR = dict(lenDicCC.items() + lenDicRR.items())
    dataListCR = alignerRobot.extractMumData(folderName, header + "Out")
    dataListCR = abunHouseKeeper.filterData(dataListCR, lenDicCR)
            
    numberOfNodes = len(lenDicCR) 
    G = graphLib.seqGraph(numberOfNodes)
    N1, N2 = len(lenDicCC), len(lenDicRR)
    print "N1, N2, numberOfNodes: ", N1, N2, numberOfNodes
    
    '''
    e.g. of dataListCC[0], dataListRR[0], dataListCR[0]
    
    [1, 520, 2913194, 2913716, 520, 523, 99.05, 'Contig0_d', 'Contig2_d']
    [1, 1383, 1253, 2603, 1383, 1351, 82.39, 'Read0_d', 'Read1705_p']
    [1, 718, 4334, 5074, 718, 741, 91.91, 'Contig0_d', 'Read1018_d']
    
    '''
    
    addDataToList(dataListCC, G, 0, 0, 'C', 'C')
    
    addDataToList(dataListRR, G, N1, N1, 'R', 'R')
    
    addDataToList(dataListCR, G, 0, N1, 'C', 'R')

    Gnew = formExtraEdges(folderName,optTypeFileHeader, contigFilename, G, N1)
    
    Gnew.saveToFile(folderName, graphName)
    
    print "len(Gnew.graphNodesList)", len(Gnew .graphNodesList)
    
    
def formExtraEdges(folderName = "/home/kakitfive/kkdata2/MetaFinisherSC/dataFolderBackup/",optTypeFileHeader="phaseString", contigFilename="improved3", G = [], N1 = 0):

    dataList = alignerRobot.extractMumData(folderName, optTypeFileHeader + "CR" + "Out")
    dataList.sort(key = itemgetter(-2))
    lenDic =  IORobot.obtainLength(folderName, contigFilename + "_Double.fasta")

    count = 0
    tmpItem = []
    embedContig2ReadDic, read2EmbedContigDic = {}, {}

    for key, items in groupby(dataList, itemgetter(-2)):
        isEmbedded = False
        for eachitem in items:
            #print eachitem
            if eachitem[4] > lenDic[key]-300 :
                    isEmbedded = True
                    tmpItem = eachitem

        if isEmbedded:
            count = count + 1
            readName = tmpItem[-1]
            embedContig2ReadDic[key] = readName
            read2EmbedContigDic[readName] = key


    print "len(embedContig2ReadDic)", len(embedContig2ReadDic)

    #assert(False)

    for contigName in embedContig2ReadDic:
        readName = embedContig2ReadDic[contigName]
    
        readIndex, contigIndex = abunHouseKeeper.parseEdgeNameToID(readName, 'R'), abunHouseKeeper.parseEdgeNameToID(contigName, 'C')

        for eachprev in G.graphNodesList[readIndex].listOfPrevNodes:
            idNode, wt = eachprev[0], eachprev[1]
            if idNode < N1:
                G.insertEdge(idNode, contigIndex, wt)

        for eachnext in G.graphNodesList[readIndex].listOfNextNodes:
            idNode, wt = eachnext[0], eachnext[1]
            if idNode < N1:
                G.insertEdge(contigIndex,idNode,  wt)


    return G


#formExtraEdges()
















