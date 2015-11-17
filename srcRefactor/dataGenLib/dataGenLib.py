# ## Generate test cases

# ## Ideal test case : two genomes with a shared repeat , with different abundances

# ## Goal : separate them 


# Have a commandline tool to generate genomes/reads/etc
import time
import argparse
import random 
from itertools import groupby
from operator import itemgetter

    
def writeListToFile(folderName, filename, seqList):
    
    f = open(folderName + filename , 'w')
    
    for i in range(len(seqList)):
        f.write(">Segkk" + str(i) + "\n")
        f.write(seqList[i] + "\n")
    
    f.close()

def randomStringGen(G):
    alphTB = ['A', 'C' , 'G', 'T']
    genome = ""
    
    for i in range(G):
        tmp = random.randint(0, 3)
        genome = genome + alphTB[tmp]
    return genome

class locList(object):
    def __init__(self):
        self.repeatInfoList = []
    def addRepeatInfo(self, repeatInfo):
        self.repeatInfoList.append(repeatInfo)


        
class repeatInfo(object):
    def __init__(self):
        self.repeatCopyList = []
        
    def addRepeatCopy(self, genomeIndex, startLoc, endLoc):
        self.repeatCopyList.append(repeatCopy(genomeIndex, startLoc, endLoc))
    
    def firstCopy(self):
        if len(self.repeatCopyList) > 0:
            return self.repeatCopyList[0]
        else:
            return None
     
        
class repeatCopy(object):
    def __init__(self, genomeIndex, startLoc, endLoc):
        self.genomeIndex = genomeIndex
        self.startLoc =  startLoc
        self.endLoc = endLoc
    
    
def insertRepeat(genomeList, locListTmp):
    
    # format : locList = [ [repeat1Information], [repeat2Information], ...    ] = [ [ [0,10,30], [1,20,40]  ] , [] , [] ,... [] ]
    # repeat1Information = [ [genomeIndex1, startLoc1, endLoc1], [] , ..., [] ] 
    
    for eachRepeatInfo in locListTmp.repeatInfoList:
        head = eachRepeatInfo.firstCopy()
        tmpSeg = genomeList[head.genomeIndex][head.startLoc:head.endLoc]
        for eachcopy in eachRepeatInfo.repeatCopyList:
            index = eachcopy.genomeIndex
            start, end = eachcopy.startLoc, eachcopy.endLoc
            genomeList[index] = genomeList[index][0:start] + tmpSeg + genomeList[index][end:]
        
    return genomeList

        
def genomeGenerate(m):
    G = 5 * pow(10, 6)
    lrepeat = pow(10, 4)
    
    startLoc = G / 2
    endLoc = startLoc + lrepeat
    
    locListTmp = locList()
    
    repeatInfoTmp = repeatInfo()
    
    for i in range(m):
        repeatInfoTmp.addRepeatCopy(i, startLoc, endLoc)
        
    locListTmp.addRepeatInfo(repeatInfoTmp)
    
    
    seqList = []
    
    for i in range(m):
        seqList.append(randomStringGen(G))
    
    
    seqList = insertRepeat(seqList, locListTmp)
    
    
    return seqList


def readFromFasta(filename):
    genList = []
    f = open(filename, 'r')
    tmp = f.readline().rstrip()
    tmpStr = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpStr) > 0:
                genList.append(tmpStr)
                tmpStr = ""
        else:
            tmpStr = tmpStr + tmp
        
        tmp = f.readline().rstrip()
    
    if len(tmpStr) > 0:
        genList.append(tmpStr)
        tmpStr = ""

    f.close()
    return genList


def contigBreakDown(contigList, breakPoints):
    ### Format : breakPoints = [ [0, 1000], [0, 2000] , [1, 299] ]
    ### [ [contigIndex, breakLoc] ]
    
    newSegList = []
    
    breakPoints.sort()
    
    numContig = len(contigList) 
    insertedList = [False for i in range(numContig)]
    
    for key, items in groupby(breakPoints, itemgetter(0)):
        startIndex = 0
        endIndex = len(contigList[0])
        contigIndex = key
        
        for eachitem in items:
            endIndex = eachitem[1]
            newSegList.append(contigList[contigIndex][startIndex:endIndex])
            startIndex = endIndex
        
        newSegList.append(contigList[contigIndex][startIndex:])
        
        insertedList[key] = True
    
    
    for i in range(numContig):
        if insertedList[i] == False:
            newSegList.append(contigList[i])
    
        
    return newSegList

    

def contigsGenerate(folderName):
    
    genList = readFromFasta(folderName + "reference.fasta")
    breakPoints = []
    
    for i in range(len(genList)):
        G = len(genList[i])
        breakPoints.append([i, G/2])
        
    segList = contigBreakDown(genList, breakPoints)
    
    return segList
    

def createANoisyRead(L, p, genome):
    G = len(genome)
    i = random.randint(0, G - L)
    segExtracted = genome[i:i + L]
    segReturn = ""
    
    
    for j in range(L):
        rand1 = random.random()
        rand2 = random.random()
        
        if rand1 >= p:
            segReturn = segReturn + segExtracted[j]
        
        if 3 * p / 4 <= rand2 < p :
            segReturn = segReturn + "A"
        elif 2 * p / 4 <= rand2 < 3 * p / 4 :
            segReturn = segReturn + "C"
        elif p / 4 <= rand2 < 2 * p / 4 :
            segReturn = segReturn + "G"
        elif rand2 < p / 4 :
            segReturn = segReturn + "T"
        
    return segReturn
    
    
    
def noisyReadsGenerate(folderName):
    print ""
    # # parameters needed : abundances, noise rate, length, number, [lambda, N, L, p]'
    c = 50
    p = 0.01
    L = 5000
    abunList = [0.4, 0.6]
    
    genList = readFromFasta(folderName + "reference.fasta")
    
    G = len(genList[0])
    N = int(G * c / L)
    
    NList = []
    segList = []
    
    for item in abunList:
        NList.append(int(item * N))    
    
    for i in range(len(genList)):
        for j in range(NList[i]):
            segList.append(createANoisyRead(L, p, genList[i]))
    
    
    return segList
    
    
def actionOnG(folderName="", numRef=""):
    print "Generating genome"
    seqList = genomeGenerate(int(numRef))
    print "Writing genome to file"
    writeListToFile(folderName , 'reference.fasta', seqList)
    

def actionOnC(folderName=""):
    print "Generating contigs"
    seqList = contigsGenerate(folderName) 
    print "Writing to file"
    writeListToFile(folderName , 'contigs.fasta', seqList)
    

def actionOnN(folderName=""):
    print "Generating noisy reads"
    seqList = noisyReadsGenerate(folderName)
    print "Writing to file"
    writeListToFile(folderName , 'raw_reads.fasta', seqList)
    
    
    
    
'''
t0 = time.time()
parser = argparse.ArgumentParser(description='Data generator')

parser.add_argument('-g', '--genomeGen', help='', required=False)
parser.add_argument('-c', '--contigsGen', help='', required=False)
parser.add_argument('-n', '--noisyReads', help='', required=False)
parser.add_argument('-p', '--folderName', help='', required=False)

args = vars(parser.parse_args())



if args['genomeGen'] != None:
    actionOnG(args['folderName'], args['genomeGen'])

if args['contigsGen'] != None:
    actionOnC(args['folderName'])
    

if args['noisyReads'] != None:
    actionOnN(args['folderName'])
    

print "args", args
print  "Time", time.time() - t0
'''


