import merger 
from ..repeatPhaserLib.finisherSCCoreLib import houseKeeper
from ..repeatPhaserLib.finisherSCCoreLib import nonRedundantResolver

import time
import argparse
import os
import json

t0 = time.time()

parser = argparse.ArgumentParser(description='mFixer')
parser.add_argument('folderName')
parser.add_argument('mummerLink')
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-t', '--LCReads', help= 'Type of reads aligned to the contigs formed by long reads', required=False)
parser.add_argument('-op', '--option', help= 'Parameters to pass in', required=False)
parser.add_argument('-cl', '--cleanInput', help= 'Parameters to clean inputs', required=False)

args = vars(parser.parse_args())

if args['parallel'] != None:
    houseKeeper.globalParallel = int(args['parallel'])
else:
    houseKeeper.globalParallel = 20

if args['LCReads'] == None:
    merger.mergerGlobalLCReads = "LR"
elif args['LCReads'] == "LR":
    merger.mergerGlobalLCReads = "LR"
elif args['LCReads'] == "SR":
    merger.mergerGlobalLCReads = "SR"

pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'] , args['mummerLink'], False)

if args['cleanInput'] == "True":
    fileList = ["SC.fasta", "LC.fasta", "SR.fasta", "LR.fasta"]

    for eachitem in fileList:
            os.system("sed -e 's/|//g' " + newFolderName + eachitem+"  > " + newFolderName + "tmp.fasta")
            os.system("cp " + newFolderName + "tmp.fasta " + newFolderName + eachitem )


if args['option'] != None:
    settingDataCombo = args['option'].split()
    settingDic = {}

    for eachitem in settingDataCombo:
        tmp = eachitem.split('=')
        settingDic[tmp[0]] = tmp[1]

    canLoad = merger.mergerGlobalFixerRobot.loadData(settingDic)
else:
    canLoad = True    


if canLoad:
    settingDic = merger.mergerGlobalFixerRobot.__dict__
    with open(newFolderName + "optionMFixer.json", 'w') as f:
        json.dump(settingDic, f)


if canLoad == True:
    merger.mainFlow(newFolderName, newMummerLink)
else: 
    print "Sorry. The above folders or files are missing or options are not correct. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"


    
print  "Time", time.time() - t0

