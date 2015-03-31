import merger 
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver

import time
import argparse
import os

t0 = time.time()

parser = argparse.ArgumentParser(description='mFixer')
parser.add_argument('folderName')
parser.add_argument('mummerLink')
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-t', '--LCReads', help= 'Type of reads aligned to the contigs formed by long reads', required=False)



args = vars(parser.parse_args())


if args['parallel'] != None:
    houseKeeper.globalParallel = int(args['parallel'])
else:
    houseKeeper.globalParallel = 1


if args['LCReads'] == None:
    merger.mergerGlobalLCReads = "SR"
elif args['LCReads'] == "LR":
    merger.mergerGlobalLCReads = "LR"
elif args['LCReads'] == "SR":
    merger.mergerGlobalLCReads = "SR"

pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'] , args['mummerLink'], False)

fileList = ["SC.fasta", "LC.fasta", "SR.fasta", "LR.fasta"]

for eachitem in fileList:
        os.system("sed -e 's/|//g' " + newFolderName + eachitem+"  > " + newFolderName + "tmp.fasta")
        os.system("cp " + newFolderName + "tmp.fasta " + newFolderName + eachitem )


if merger.mergerGlobalLCReads == "SR":
    merger.mergeContigs(newFolderName , newMummerLink)
    command = "cp "+ newFolderName + "LR.fasta "+ newFolderName + "raw_reads.fasta"
    os.system(command)
    print "Command: ",  command 
    
elif merger.mergerGlobalLCReads == "LR":
    
    command = "cp "+ newFolderName + "LR.fasta " + newFolderName + "SR.fasta "
    os.system(command)
    print "Command: ",  command 
    
    merger.onlyLRMiassemblyFix(newFolderName, newMummerLink)
    
    command = "cp "+ newFolderName + "LC_n.fasta "+ newFolderName + "contigs.fasta"
    os.system(command)
    print "Command: ",  command 
    
    command = "cp "+ newFolderName + "LR.fasta "+ newFolderName + "raw_reads.fasta"
    os.system(command)
    print "Command: ",  command 
    
    nonRedundantResolver.removeEmbedded(newFolderName , newMummerLink)

    
    
    
    
    
    
print  "Time", time.time() - t0

