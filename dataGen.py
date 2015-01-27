# ## Generate test cases

# ## Ideal test case : two genomes with a shared repeat , with different abundances

# ## Goal : separate them 


# Have a commandline tool to generate genomes/reads/etc
import time
import argparse
import random 

G = 5 * pow(10, 6)
lrepeat = pow(10, 4)
    
    
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

def genomeGenerate(m):
    
    startLoc = G / 2
    endLoc = startLoc + lrepeat
    seqList = []
    
    for i in range(m):
        seqList.append(randomStringGen(G))
    
    repeatSeg = seqList[0][startLoc:endLoc]
     
    for i in range(m):
        tmp = seqList[i][0:startLoc] + repeatSeg + seqList[i][endLoc:]
        seqList[i] = tmp
    
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
    
def contigsGenerate(folderName):
    print ""
    segList = []
    genList = readFromFasta(folderName + "reference.fasta")
    
    for eachgen in genList:
        segList.append(eachgen[0:G / 2])
        segList.append(eachgen[G / 2:])

    return segList
    

def createANoisyRead(L, p, genome):
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
    
    totalLen = 0
    for eachitem in genList:
        totalLen = totalLen + len(eachitem)
    
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


