'''
Specification : 

	a) python setPara.py --outMFixer mFixerAttr.txt --setting runa=True ...
	b) python setPara.py --outASplitter aSplitterAttr.txt --setting runa=True...

'''

import time
import argparse
import os
import json

from repeatPhaserLib import abunHouseKeeper
from repeatPhaserLib.finisherSCCoreLib import houseKeeper

t0 = time.time()

parser = argparse.ArgumentParser(description='setPara')

parser.add_argument('folderName')
parser.add_argument('mummerLink')

parser.add_argument('--outMFixer' , required=False)
parser.add_argument('--outASplitter' , required=False)



args = vars(parser.parse_args())

pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'], args['mummerLink'])



if args['outASplitter'] != None:
	settingDataCombo = args['setting'].split()
	settingDic = {}

	for eachitem in settingDataCombo:
		tmp = eachitem.split('=')
		settingDic[tmp[0]] = tmp[1]

	with open(newFolderName + args['outASplitter'], 'w') as f:
		json.dump(settingDic, f)


#abunGlobalSplitParameterRobot