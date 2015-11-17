import time
import datetime
import argparse
import os
import json
import re


from repeatPhaserLib import abunHouseKeeper
from repeatPhaserLib.finisherSCCoreLib import houseKeeper

from misassemblyFixerLib import merger
	
def runEvaluation(quastPath, folderName,paraMFixerFileName,paraASplitterFileName, outputFilename):
	command = "python "+quastPath+" " + folderName +"abun.fasta -o "  + folderName + "   -R  " + folderName + "reference.fasta"
	os.system(command)
	
	f = open(outputFilename , 'a')
	tmpString  = ""
	st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d~%H:%M:%S')
	
	tmpString = tmpString + str(st) + ","

	with open(folderName+paraMFixerFileName) as fdicMFixer:
		myDicMFixer = json.load(fdicMFixer)

	for eachitem in myDicMFixer:
		tmpString = tmpString + str(myDicMFixer[eachitem]) + ","
    

	with open(folderName+paraASplitterFileName) as fdicASplitter:
		myDicASplitter = json.load(fdicASplitter)


	for eachitem in myDicASplitter:
		tmpString = tmpString + str(myDicASplitter[eachitem]) + ","
	
	quastReport = folderName + "report.txt"
	
	regex = re.compile("misassemblies*")

	misassembly = str(1997)
	ncontig = str(1997) 
	
	with open(quastReport) as f2:
	    for line in f2:
	        result = re.match("# misassemblies (.*?) .*", line)
	        if result != None: 
		        ans = result.group(0).split()
		        misassembly = ans[2]
	        
	        result = re.match("# contigs \(>= 0 bp\) (.*?) .*", line)
	        if result != None: 
		        ans = result.group(0).split()
		        ncontig=ans[-1]

	
	tmpString = tmpString  + misassembly + "," + ncontig+ "\n"

	f.write(tmpString)

	f.close()


def loggHeaders(folderName, outputFilename):

	myDicMFixer = merger.mergerGlobalFixerRobot.__dict__
	myDicASplitter =  abunHouseKeeper.abunGlobalSplitParameterRobot.__dict__
	
	f = open(outputFilename , 'a')
	tmpString  = ""
	st = "time"
	
	tmpString = tmpString + str(st) + ","

	for eachitem in myDicMFixer:
		tmpString = tmpString + eachitem + ","

	for eachitem in myDicASplitter:
		tmpString = tmpString + eachitem + ","


	tmpString = tmpString  + "misassembly" + "," + "ncontig"+ "\n"

	f.write(tmpString)
	f.close()


parser = argparse.ArgumentParser(description='evaluation')

parser.add_argument('--option', required = True)
parser.add_argument('--quastPath' , required=False)
parser.add_argument('--folderName' , required=False)
parser.add_argument('--paraMFixerFileName' , required=False)
parser.add_argument('--paraASplitterFileName' , required=False)
parser.add_argument('--outputFilename' , required=False)


args = vars(parser.parse_args())
quastPath, folderName,paraMFixerFileName,paraASplitterFileName, outputFilename = args['quastPath'], args['folderName'], args['paraMFixerFileName'],args['paraASplitterFileName'], args['outputFilename']

if args['option'] == 'header' :
	loggHeaders(folderName, outputFilename)
elif args['option'] == 'evaluate':
	runEvaluation(quastPath, folderName,paraMFixerFileName,paraASplitterFileName, outputFilename)







