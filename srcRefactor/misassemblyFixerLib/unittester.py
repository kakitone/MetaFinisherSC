
import unittest
import os
import merger
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper
from ..repeatPhaserLib import repeatFlankingDefiner
from ..repeatPhaserLib.finisherSCCoreLib import IORobot


class IsOddTests(unittest.TestCase):
    
    def setUp(self):
        print "Init : Copying files : "
        self.testingFolder = "unitTestFolder"
        self.mummerPath = "/Users/kakitlam/Desktop/experimentBench/MUMmer3.23/"
        self.listOfFiles = ["raw_reads.fasta", "contigs.fasta"]
        os.system("rm -rf "+ self.testingFolder)
        
        
    def testLC_SR(self):
        print "Testing LC_SR "
        
        if True:
            self.runningLCSRTestSet("/Users/kakitlam/Desktop/metaFinisherSC/LC_SR/")
        
            matchingContigFile = "LC_n.fasta"
            ctexpected = 4
            lenDic = IORobot.obtainLength(self.testingFolder +"/",matchingContigFile  )
            assert(len(lenDic) == ctexpected)
        
            os.system("rm -rf "+ self.testingFolder)
    
    def testSC_LR(self):
        print "Testing SC_LR "
        if True:
            self.runningSCLRTestSet("/Users/kakitlam/Desktop/metaFinisherSC/SC_LR/")
        
            matchingContigFile = "SC_n.fasta"
            ctexpected = 4
            
            lenDic = IORobot.obtainLength(self.testingFolder+"/",  matchingContigFile)
            assert(len(lenDic) == ctexpected)
            os.system("rm -rf "+ self.testingFolder)
    
    def testMerger(self):
        print "Testing SC_LR and LC_SR "
        
        if True:
            self.runningMergerTestSet()
            matchingContigFile = "contigs.fasta"
            ctexpected = 8
            lenDic = IORobot.obtainLength(self.testingFolder+"/",  matchingContigFile)
            print lenDic
            print "len(lenDic)", len(lenDic)
            assert(len(lenDic) == ctexpected)
            
            os.system("rm -rf "+ self.testingFolder)
        
    
    def testIntegrationTestLC_SR(self):
        print "Integration Test LC_SR, FinisherSC and ASplitter"
        # Current focus : 
        if True:
        
            self.runningLCSRTestSet("/Users/kakitlam/Desktop/metaFinisherSC/LC_SR/")
        
            
            os.system("cp "+ self.testingFolder + "/LC_n.fasta " + self.testingFolder + \
                       "/contigs.fasta")
            
            
            os.system("cp "+ self.testingFolder + "/SR.fasta " + self.testingFolder + \
                      "/raw_reads.fasta")
            
            os.system("python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC " +\
                      "-par 4 "+ self.testingFolder + " "+ self.mummerPath)
            
            
            os.system("python -m srcRefactor.repeatPhaserLib.aSplitter -par 4 " + \
                      self.testingFolder + " " + self.mummerPath )
            
            matchingContigFile = "abun.fasta"
            ctexpected = 3
            
            lenDic = IORobot.obtainLength(self.testingFolder+"/",  matchingContigFile)
            assert(len(lenDic) == ctexpected)
            
            os.system("rm -rf "+ self.testingFolder)
        
        
    def testIntegrationTestSC_LR(self):
        print "Integration Test SC_LR, FinisherSC"    
        
        if True:
            self.runningSCLRTestSet("/Users/kakitlam/Desktop/metaFinisherSC/SC_LR/")
            
                
            os.system("cp "+ self.testingFolder + "/SC_n.fasta " + self.testingFolder + \
                       "/contigs.fasta")
            
            
            os.system("cp "+ self.testingFolder + "/LR.fasta " + self.testingFolder + \
                      "/raw_reads.fasta")
            
            
            os.system("python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC " +\
                      "-par 4 -th True "+ self.testingFolder + " "+ self.mummerPath)
            

            matchingContigFile = "improved2.fasta"
            ctexpected = 1
            
            lenDic = IORobot.obtainLength(self.testingFolder+"/",  matchingContigFile)
            assert(len(lenDic) == ctexpected)
            
        
            os.system("rm -rf "+ self.testingFolder)
        
            
    def testFullIntegrationTest(self):
        print "Full Integration Test on SC, LC, LR, SR ; " +\
              "using Merger, FinisherSC, ASplitter."
        
        if True:
            self.runningMergerTestSet()
        
            os.system("cp "+ self.testingFolder + "/LR.fasta " + self.testingFolder + \
                "/raw_reads.fasta")
      
            os.system("python -m srcRefactor.repeatPhaserLib.finisherSCCoreLib.finisherSC " +\
                      "-par 4 "+ self.testingFolder + " "+ self.mummerPath)
            
    
            os.system("python -m srcRefactor.repeatPhaserLib.aSplitter -par 4 " + \
                      self.testingFolder + " " + self.mummerPath )
            
            matchingContigFile = "abun.fasta"
            ctexpected = 4
            
            lenDic = IORobot.obtainLength(self.testingFolder+"/",  matchingContigFile)
            assert(len(lenDic) == ctexpected)
            
            os.system("rm -rf "+ self.testingFolder)
    
    def runningMergerTestSet(self):
        self.listOfFiles = ["/Users/kakitlam/Desktop/metaFinisherSC/LC_SR/LC.fasta", \
                            "/Users/kakitlam/Desktop/metaFinisherSC/LC_SR/SR.fasta", \
                            "/Users/kakitlam/Desktop/metaFinisherSC/SC_LR/SC.fasta", \
                            "/Users/kakitlam/Desktop/metaFinisherSC/SC_LR/LR.fasta"]
        
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp " + eachitem + " " +self.testingFolder)
        
        os.system("python -m srcRefactor.misassemblyFixerLib.mFixer -par 4 "+ self.testingFolder + "/ "+ self.mummerPath )
        
        
        
    def runningLCSRTestSet(self ,myFolderName):
        print "Integration test on Merger:  " + myFolderName
        self.listOfFiles = ["LC.fasta", "SR.fasta"]
        self.sourceFolder = myFolderName
        
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp "+ self.sourceFolder + eachitem + " " +self.testingFolder)
        
        houseKeeper.globalParallel = 4
        merger.fixLCMisassembly(self.testingFolder+"/", self.mummerPath)
        
        #os.system("python -m srcRefactor.misassemblyFixerLib.mFixer -par 4 -t LR "+ self.testingFolder + "/ "+ self.mummerPath )
        
    def runningSCLRTestSet(self, myFolderName):
        print "Integration test on Merger:  " + myFolderName
        self.listOfFiles = ["SC.fasta", "LR.fasta", "LC_n.fasta"]
        self.sourceFolder = myFolderName
        
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp "+ self.sourceFolder + eachitem + " " +self.testingFolder)
        
        houseKeeper.globalParallel = 4
        merger.fixSCMisassembly(self.testingFolder+"/", self.mummerPath)
        
        
    def tearDown(self):
        print "Teardown : Removing used files "
        
        

def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
