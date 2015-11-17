import json
import numpy as np
import os
from itertools import groupby
from operator import itemgetter
import copy
import sys
from random import shuffle
from multiprocessing import Pool

from finisherSCCoreLib import alignerRobot
from finisherSCCoreLib import graphLib
from finisherSCCoreLib import IORobot
from finisherSCCoreLib import houseKeeper
from finisherSCCoreLib import nonRedundantResolver 

import associatedReadFinder
import readContigGraphFormer
import repeatFinder
import abunHouseKeeper
import abunGraphLib
import emalgo

def mainFlow(folderName, mummerLink):

    contigFilename = "improved3"
    readsetFilename = "phasingSeedName"
    optTypeFileHeader = "phaseString"
    contigReadGraph = "phaseStringGraph1"
    repeatFilename = "phaseRepeat.txt"
    repeatSpec = "repeatSpecification.txt"
    optionToRun = "xphase"

    if False:
		associatedReadFinder.getAllAssociatedReads(folderName, mummerLink,readsetFilename)
   		readContigGraphFormer.formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph )
    	repeatFinder.identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun )
    
    if False:
	    myCountDic = generateAbundanceGraph(folderName, mummerLink, contigFilename)
    
	if False:
    	splitter(folderName, mummerLink, contigReadGraph, contigFilename,readsetFilename )

        
	os.system("cp selected_raw.part-* "+ folderName )
	os.system("rm selected_raw.part-*")

#mainFlow()