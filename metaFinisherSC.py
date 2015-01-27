'''
import os
from itertools import groupby
from operator import itemgetter
import sys
'''

import time
import argparse
from newPhasingMain import newPhasing
from abundanceEstimation import abun
import json
import os

def splitter(folderName, mummerLink):

    '''
    Input : repeatSpecification.txt , myCountDic.json, improved3.fasta, raw_reads.fasta
    Output : abunsplit.fasta
    
    Algorithm : 
    
    1. Load data from various sources [various json files]
    
    2. For each repeat interior:
        a) identify the abundances associated with in/out contigs
        b) perform a split and record the split
    
    3. Use split results to generate contigs [may already exist in newPhasing.py ] 
        a) use a graph to capture the split results 
        b) use reads to fill in any gaps 
        c) read out the contigs 
    
    '''
    
    
    with open(folderName + 'myCountDic.json') as f:
        myCountDic = json.load(f)
        
    newPhasing.abunSplit(folderName, mummerLink, myCountDic)    
    
    
def mainFlow(folderName, mummerLink):
    print "Hello world"
    
    contigFilename = "improved3"
    readsetFilename = "phasingSeedName"
    optTypeFileHeader = "phaseString"
    contigReadGraph = "phaseStringGraph1"
    repeatFilename = "phaseRepeat.txt"
    repeatSpec = "repeatSpecification.txt"
    optionToRun = "xphase"
    
    os.system("python finisherSC.py "+ folderName +" " + mummerLink)
    
    if True:
        newPhasing.getAllAssociatedReads(folderName, mummerLink,readsetFilename)
        newPhasing.formReadContigStringGraph(folderName, mummerLink,contigFilename, readsetFilename, optTypeFileHeader , contigReadGraph )
        newPhasing.identifyRepeat(folderName, mummerLink,contigFilename,contigReadGraph, repeatFilename, optionToRun )
        newPhasing.defineRepeatAndFlanking(folderName, mummerLink,contigFilename,contigReadGraph,repeatFilename,repeatSpec )
        
    if True:
        myCountDic = abun.generateAbundanceGraph(folderName, mummerLink)
        
    if True :
        splitter(folderName, mummerLink)
        
t0 = time.time()

parser = argparse.ArgumentParser(description='MetaFinisherSC')
parser.add_argument('folderName')
parser.add_argument('mummerLink')

args = vars(parser.parse_args())
print "args", args
mainFlow(args['folderName'] , args['mummerLink'])
print  "Time", time.time() - t0


