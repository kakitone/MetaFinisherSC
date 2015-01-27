import matplotlib.pyplot as plt
import commonLib
import random
import time
import os
from itertools import groupby
from operator import itemgetter
import copy
import json

class seqGraphNodeWt(commonLib.seqGraphNode):
    def __init__(self, nodeIndex):
        commonLib.seqGraphNode.__init__(self, nodeIndex)
        self.nodeWt = 0 
        
class seqGraphWt(commonLib.seqGraph):
    def __init__(self, numberOfNodes):
        self.graphNodesList = [seqGraphNodeWt(i) for i in range(numberOfNodes)] 
    def formReportName(self, index):
        return str(index) + "_" +str(self.graphNodesList[index].nodeWt)


'''
Input : string_graph_3, improved3.fasta, raw_reads.fasta
Output : string_graph_4 with weights [need a data structure to store the weight on node]

Algorithm : 
1. Find your favorite mappers to map read back
    a. MUMmer, Bowtie, bbmap, any that works V 
    b. And then write a short parser to parse the results V 
2. Calculate count on the abundances 
    a. Aggregate by taking average [put weights on bin along contigs]
    b. Inheritance and a subclass 
3. Find your favorite graphical tool to display 
    a. Use a javascript library [halfviz should just work ! put weight on edge ]

'''

def evaluateCoverage(dataList, lenDic, readLenDic, folderName,mummerLink, continueFilter):
    
    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = 0
            
    dataList.sort(key = itemgetter(-1)) 
    
    ctkk, ctbase = 0, 0
    toAddBackDic = copy.deepcopy(readLenDic)
    
    for key, items in groupby(dataList, itemgetter(-1)):
        maxMatch = -1
        bestname = ""
        
        for eachitem in items:
            ct = eachitem[6]/100.0 * eachitem[4]
            if ct > maxMatch:
                maxMatch = ct 
                bestname = eachitem[-2]
        myCountDic[bestname] += readLenDic[key] 
        
        ctkk = ctkk + 1 
        ctbase = ctbase + readLenDic[key]
        toAddBackDic[key] = -1
    
    cttot = 0
    for eachitem in readLenDic:
        cttot = cttot + readLenDic[eachitem]
        
    print "Missed coverage  ", (cttot - ctbase)/(4.7*pow(10, 6))
    print "percentage miss read", (len(readLenDic) - ctkk)/(1.0*len(readLenDic)) 
    
    toAddReadList = []
    for eachitem in toAddBackDic:
        if toAddBackDic[eachitem] >= 0 :
            toAddReadList.append(eachitem)
    
    
    if continueFilter:
        commonLib.putListToFileO(folderName, "raw_reads.fasta" , "selected_raw", toAddReadList)
        
        command = mummerLink + "nucmer -l 10 --maxmatch -p " + folderName + "kkout " + folderName + "improved3.fasta " + folderName + "selected_raw.fasta"
        os.system(command)
        
        command = mummerLink + "show-coords -r " + folderName + "kkout.delta > " + folderName + "abunMissOut"
        os.system(command)
        
    for i in range(len(myCountDic)):
        eachitem = "Segkk"+str(i)
        print eachitem , myCountDic[eachitem]/(1.0*lenDic[eachitem])
        myCountDic[eachitem] = myCountDic[eachitem]/(1.0*lenDic[eachitem])
        
    return myCountDic
    
def generateAbundanceGraph(folderName, mummerLink):
    
    
    print "generateAbundanceGraph"
    
    '''
    1. Find your favorite mappers to map read back
        a. MUMmer, Bowtie, bbmap, any that works V 
        b. And then write a short parser to parse the results V 
    '''
    numberOfFiles = 20
    
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
    
        if True:
            command = mummerLink + "nucmer --maxmatch --nosimplify -p " + folderName + "out " + folderName + "improved3.fasta "+folderName+"raw_reads.part-" + indexOfMum + ".fasta"
            os.system(command)
        
            command = mummerLink + "show-coords -r " + folderName + "out.delta > " + folderName + "fromMumAbun" + indexOfMum
            os.system(command)
        
    dataList = []
    
    for i in range(1, 1+numberOfFiles):
        if i < 10:
            indexOfMum = "0" + str(i)
        else:
            indexOfMum = str(i)
        dataList = dataList+ commonLib.extractMumData(folderName, "fromMumAbun"+ str(indexOfMum))

    '''
    2. Calculate count on the abundances 
        a. Aggregate by taking average [put weights on bin along contigs]
        b. Inheritance and a subclass 
    '''
         
    lenDic = commonLib.obtainLength(folderName, "improved3.fasta")
    readLenDic = commonLib.obtainLength(folderName , "raw_reads.fasta")
    

    myCountDic = {}
    for eachitem in lenDic:
        myCountDic[eachitem] = [0 for i in range(lenDic[eachitem])]

    thres = 30
    lenSum = 0
    extraDataList= []
    
    
    print "len(dataList)", len(dataList)
    myCountDic =  evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,  True)
    extraDataList = commonLib.extractMumData(folderName, "abunMissOut" )
    dataList = dataList + extraDataList
    myCountDic = evaluateCoverage(dataList, lenDic, readLenDic, folderName, mummerLink,False)
    
    with open(folderName + 'myCountDic.json', 'w') as f:
        json.dump(myCountDic, f)

    
    return myCountDic
    
'''
folderName, mummerLink = "/Users/kakitlam/Desktop/Research/Dec15-2014/src/EcoliTestRun/", "/Users/kakitlam/Desktop/Research/Dec15-2014/src/MUMmer3.23/"
t0 = time.time()
generateAbundanceGraph(folderName, mummerLink)
print "Time : seconds ", time.time() - t0, time.time()-t0 -30
'''