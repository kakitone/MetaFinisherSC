
import unittest
import os
from srcRefactor.repeatPhaserLib.finisherSCCoreLib import IORobot


class IsOddTests(unittest.TestCase):
    
    def setUp(self):
        print "Init : Copying files : "
        #os.system("python ~/experimentBench/test.py")
        
        self.testingFolder = "unitTestFolder"
        self.mummerPath = "/tmp/MUMmer3.23/"
        self.listOfFiles = ["LC.fasta", "LR.fasta"]
        os.system("rm -rf "+ self.testingFolder)
        
    # Remember to have test* at the beginning

        
    def testMetaFinisherSCSyntheticTest(self):
        print "Integration test on MetaFinisherSC:  "
        self.runningTestSet("/tmp/syntheticrun/" , 2)

   
    def runningTestSet(self ,myFolderName, ctexpected):
        print "Integration test on MetaFinisherSC:  " + myFolderName
        self.sourceFolder = myFolderName
        os.system("mkdir " + self.testingFolder)
        
        for eachitem in self.listOfFiles:
            os.system("cp "+ self.sourceFolder + eachitem + " " +self.testingFolder)
        
        os.system("python -m srcRefactor.misassemblyFixerLib.mFixer -par 4 "+ self.testingFolder + " "+ self.mummerPath)
        os.system("python -m srcRefactor.repeatPhaserLib.aSplitter -par 4 "+ self.testingFolder + " "+ self.mummerPath)
        #os.system("python -m srcRefactor.repeatPhaserLib.aSplitter --pickup split "+ self.testingFolder + " "+ self.mummerPath)



        lenDic = IORobot.obtainLength(self.testingFolder, "/abun.fasta")
        print lenDic
        assert(len(lenDic) == ctexpected)
        #os.system("rm -rf "+ self.testingFolder)
        
    def tearDown(self):
        print "Teardown : Removing used files "
        
        

def main():
    unittest.main()
    
if __name__ == '__main__':
    main()
