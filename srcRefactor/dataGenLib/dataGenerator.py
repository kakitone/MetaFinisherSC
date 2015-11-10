import dataGenLib
import os

'''

Goal : This simple script generates data needed for testing. 

Input parameters : 
1) Reference
    a) G
    b) Number of vectors
    c) Any repeats ? 
    d) Repeat number ?
    e) More repeats ? 
    f) Repeat locations 
    g) Circular or not ? 
    
2) Contigs
    a) Specify how to extract it from the reference
        in particular the break locations 
    
3) Raw reads
    a) noise level
    b) noise type
    c) length L 
    d) coverage  on each genome on reference
    
Key outputs are : 
1) reference.fasta
2) contigs.fasta
3) raw_reads.fasta 


Dependency : 
1) It will incorporate several command line tools to do the job for you. 
2) It mainluy glues things together


'''



class dataGeneartorRobot(object):
    def __init__(self, folderName = ""):
        self.folderName =  folderName
        os.system("mkdir "+ folderName)
        self.setParameters()
    def generateReference(self):
        print
    def generateContigs(self):
        print
    def genearteReads(self):
        genList = dataGenLib.readFromFasta(self.folderName + "reference.fasta")
        
        NList = []
        segList = []
        
        for c in self.abun:
            NList.append(int(self.G*c/self.L))    
        
        for i in range(len(genList)):
            for j in range(NList[i]):
                segList.append(dataGenLib.createANoisyRead(self.L, self.p, genList[i]))
        
        dataGenLib.writeListToFile(self.folderName, "raw_reads.fasta", segList)
        
    def setParameters(self):
        print
    def runAll(self):
        print "generateReference"
        self.generateReference()
        
        print "generateContigs"
        self.generateContigs()
        
        print "genearteReads"
        self.genearteReads()
        

class twoSplit(dataGeneartorRobot):
    def generateReference(self):
        dataGenLib.actionOnG(self.folderName, 2)  
    def generateContigs(self):
        dataGenLib.actionOnC(self.folderName)
    def genearteReads(self):
        dataGenLib.actionOnN(self.folderName)
    

class LC_SR_example(dataGeneartorRobot):
        
    def generateReference(self):
        lrepeat = 5000
    
        startLoc = self.G / 2
        endLoc = startLoc + lrepeat
        
        locListTmp = dataGenLib.locList()
        
        repeatInfoTmp = dataGenLib.repeatInfo()
        
        for i in range(self.m):
            repeatInfoTmp.addRepeatCopy(i, startLoc, endLoc)
            
        locListTmp.addRepeatInfo(repeatInfoTmp)
        
        
        seqList = []
        
        for i in range(self.m):
            seqList.append(dataGenLib.randomStringGen(self.G))
        
        seqList = dataGenLib.insertRepeat(seqList, locListTmp)
        
        dataGenLib.writeListToFile(self.folderName , 'reference.fasta', seqList)
        
        
    def generateContigs(self):
        
        genList = dataGenLib.readFromFasta(self.folderName + "reference.fasta")
        breakPoints = []
    
        for i in range(len(genList)):
            G = len(genList[i])
            breakPoints.append([i, G/2])
            
        segList = dataGenLib.contigBreakDown(genList, breakPoints)
        
        for eachitem in segList:
            print len(eachitem)
            
        seqListNew = [segList[0] + segList[3], segList[2]+ segList[1]]
    
        dataGenLib.writeListToFile(self.folderName , 'contigs.fasta', seqListNew)
    
        
    def setParameters(self):
        self.G= 5*(10**6)
        self.L =100
        self.p = 0.01
        self.m = 2
        self.abun = [50, 20]
        

class SC_LR_example(dataGeneartorRobot):
    
    def generateReference(self):
        lrepeat = 1000
        numCopy = 3
        separation = 10**6
        
        locListTmp = dataGenLib.locList()
        repeatInfoTmp = dataGenLib.repeatInfo()
        
        for i in range(1, 1+ numCopy):
            repeatInfoTmp.addRepeatCopy(0, i*separation, i*separation + lrepeat)
        
        locListTmp.addRepeatInfo(repeatInfoTmp)
        
        seqList = [dataGenLib.randomStringGen(self.G)]        
        seqList = dataGenLib.insertRepeat(seqList, locListTmp)
        
        dataGenLib.writeListToFile(self.folderName , 'reference.fasta', seqList)
         
    def generateContigs(self):
        genList = dataGenLib.readFromFasta(self.folderName + "reference.fasta")
        breakPoints = []
        
        numCopy = 3
        separation = 10**6
        
        for i in range(1 , 1+numCopy):
            breakPoints.append([0, i*separation])
            
        segList = dataGenLib.contigBreakDown(genList, breakPoints)
        
        seqListNew = [segList[0] + segList[2]+segList[1]+ segList[3]]
    
        dataGenLib.writeListToFile(self.folderName , 'contigs.fasta', seqListNew)


    def setParameters(self):
        self.G= 5*(10**6)
        self.L = 5000
        self.p = 0.1
        self.m = 1
        self.abun = [50]
        

### Entrance point


# dataRobot = twoSplit("/tmp/test/")


#dataRobot = LC_SR_example("/tmp/LC_SR/")
#dataRobot.runAll()

#dataRobot = SC_LR_example("/tmp/SC_LR/")
#dataRobot.runAll()



### Entrance point


# dataRobot = twoSplit("/tmp/test/")


#dataRobot = LC_SR_example("/tmp/LC_SR/")
#dataRobot.runAll()

#dataRobot = SC_LR_example("/tmp/SC_LR/")
#dataRobot.runAll()




class abunGap_example(dataGeneartorRobot):
        
    def generateReference(self):
        lrepeat = self.repeatLen
    
        startLoc = self.G / 2
        endLoc = startLoc + lrepeat
        
        locListTmp = dataGenLib.locList()
        
        repeatInfoTmp = dataGenLib.repeatInfo()
        
        for i in range(self.m):
            repeatInfoTmp.addRepeatCopy(i, startLoc, endLoc)
            
        locListTmp.addRepeatInfo(repeatInfoTmp)
        
        
        seqList = []
        
        for i in range(self.m):
            seqList.append(dataGenLib.randomStringGen(self.G))
        
        seqList = dataGenLib.insertRepeat(seqList, locListTmp)
        
        dataGenLib.writeListToFile(self.folderName , 'reference.fasta', seqList)
        
        
    def generateContigs(self):
        
        genList = dataGenLib.readFromFasta(self.folderName + "reference.fasta")
        breakPoints = []
    
        for i in range(len(genList)):
            G = len(genList[i])
            breakPoints.append([i, G/2])
            
        segList = dataGenLib.contigBreakDown(genList, breakPoints)
        
        seqListNew = [segList[0] + segList[3], segList[2]+ segList[1] ]
    
    
        for eachitem in seqListNew:
            print len(eachitem)
        
        dataGenLib.writeListToFile(self.folderName , 'contigs.fasta', seqListNew)
        
        
    def setParameters(self):
        self.G= 5*(10**6)
        self.L = 6000
        self.p = 0.01
        self.m = 2
        self.abun = [50, 20]
        self.repeatLen = 12000
        

#dataRobot = abunGap_example("/tmp/abunGap/")
#dataRobot.runAll()

class polyGap_example(dataGeneartorRobot):
    def generateReference(self):
        lrepeat = self.repeatLen
    
        startLoc = self.G / 2
        endLoc = startLoc + lrepeat
        
        locListTmp = dataGenLib.locList()
        
        repeatInfoTmp = dataGenLib.repeatInfo()
        
        for i in range(self.m):
            repeatInfoTmp.addRepeatCopy(i, startLoc, endLoc)
            
        locListTmp.addRepeatInfo(repeatInfoTmp)
        
        
        seqList = []
        
        for i in range(self.m):
            seqList.append(dataGenLib.randomStringGen(self.G))
        
        seqList = dataGenLib.insertRepeat(seqList, locListTmp)

        newSeqList = ["", ""]

        pt1, pt2 = startLoc + lrepeat/3, startLoc + 2*lrepeat/3
        print pt1, pt2
        newSeqList[0] =  seqList[0][0:pt1] + 'A' +seqList[0][pt1+1:pt2] +'G' +seqList[0][pt2+1:] 
        newSeqList[1] =  seqList[1][0:pt1] + 'C' +seqList[1][pt1+1:pt2] +'T' +seqList[1][pt2+1:] 
        
        dataGenLib.writeListToFile(self.folderName , 'reference.fasta', newSeqList)
        
        
    def generateContigs(self):
        print "Generate contigs"
        genList = dataGenLib.readFromFasta(self.folderName + "reference.fasta")
        breakPoints = []
    
        for i in range(len(genList)):
            G = len(genList[i])
            breakPoints.append([i, G/2])
            
        segList = dataGenLib.contigBreakDown(genList, breakPoints)
        
        seqListNew = [segList[0] + segList[3], segList[2]+ segList[1] ]
    
    
        for eachitem in seqListNew:
            print len(eachitem)
        
        dataGenLib.writeListToFile(self.folderName , 'contigs.fasta', seqListNew)


    def setParameters(self):
        self.G= 5*(10**6)
        self.L = 6000
        self.p = 0.01
        self.m = 2
        self.abun = [50, 50]
        self.repeatLen = 12000
        


dataRobot = polyGap_example("/tmp/polyGap/")
dataRobot.runAll()

