import abunSplitter
import abunGraphLib

from finisherSCCoreLib import houseKeeper
from finisherSCCoreLib import alignerRobot

import time
import argparse
import abunHouseKeeper
import os 
import json

t0 = time.time()

parser = argparse.ArgumentParser(description='aSplitter')
parser.add_argument('folderName')
parser.add_argument('mummerLink')

parser.add_argument('-f', '--fast', help= 'Fast aligns contigs (input is True)', required=False)
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-l', '--large', help= 'Large number of contigs/large size of contigs (input is True)', required=False)

parser.add_argument('-rp', '--replace', help= 'Input files to aSplitter(e.g. noEmbed.fasta, improved.fasta, improved2.fasta or improved3.fasta)', required=False)
parser.add_argument('-ar', '--avoidrefine', help= 'Avoid refined abundance estimation (input is True)', required=False)
parser.add_argument('-rs', '--readsearch', help= 'Number of linking reads across a gap  (input is number of such linking reads/2)', required=False)
parser.add_argument('-rd', '--RRDisable', help= 'Whether one should disable Read to Read overlap check (input is True)', required=False)
parser.add_argument('-pk', '--pickup', help= 'where to run ASplitter, map/count/split/graph', required=False)

parser.add_argument('-em', '--runemalgo', help= 'whether to run EM for splitting repeats', required=False)
parser.add_argument('-mpi', '--runmpi', help= 'whether to run MPI', required=False)
parser.add_argument('-sc', '--segmentcount', help= 'number of segments to break', required=False)



parser.add_argument('-op', '--option', help='File of parameter list (input is opa=true opb=false)', required=False)

args = vars(parser.parse_args())



if args['runmpi'] == "True":
    houseKeeper.globalRunMPI = True
else:
    houseKeeper.globalRunMPI = False


me = 0 

if houseKeeper.globalRunMPI == True:
    from mpi4py import MPI
    from mpi4py.MPI import ANY_SOURCE
    comm = MPI.COMM_WORLD
    me = comm.Get_rank()
    nproc = comm.Get_size()

if me ==0 :
    print "args", args
    pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'], args['mummerLink'])

    if args['segmentcount'] != None:
        houseKeeper.globalParallelFileNum = int(args['segmentcount'])
    else:
        houseKeeper.globalParallelFileNum = 20

    if args['fast'] == "True":
        houseKeeper.globalFast = True
    else:
        houseKeeper.globalFast = False

    if args['parallel'] != None:
        houseKeeper.globalParallel = int(args['parallel'])
    else:
        houseKeeper.globalParallel = 20


    if args['large'] == "True":
        houseKeeper.globalLarge = True
    else:
        houseKeeper.globalLarge = False


    if args['avoidrefine'] == "False":
        abunHouseKeeper.abunGlobalAvoidrefine = False
    else:
        abunHouseKeeper.abunGlobalAvoidrefine = True


    if args['readsearch'] != None:
        abunHouseKeeper.abunGlobalReadSearchDepth = int(args['readsearch']) 
    else:
        abunHouseKeeper.abunGlobalReadSearchDepth = 0


    if args['replace'] != None : 
        if  args['replace'] == 'skip':
            print "skip copy"
        else:
            abunHouseKeeper.replaceFiles( newFolderName, args['replace']) 
    else:
        abunHouseKeeper.replaceFiles( newFolderName, "mFixed.fasta")

    if args['RRDisable'] == "False":
        abunHouseKeeper.abunGlobalRRDisable = False
    else:
        abunHouseKeeper.abunGlobalRRDisable = True


    if args['pickup'] in [ "map", "count", "split", "graph"] :
        abunHouseKeeper.abunGlobalRunPickUp = args['pickup']

    if args['runemalgo'] == "True":
        abunHouseKeeper.abunGlobalRunEM = True
    else:
        abunHouseKeeper.abunGlobalRunEM = False

    if args['option'] != None:
        settingDataCombo = args['option'].split()
        settingDic = {}

        for eachitem in settingDataCombo:
            tmp = eachitem.split('=')
            settingDic[tmp[0]] = tmp[1]

        canLoad = abunHouseKeeper.abunGlobalSplitParameterRobot.loadData(settingDic)
    else:
        canLoad = True    

    if canLoad:
        settingDic = abunHouseKeeper.abunGlobalSplitParameterRobot.__dict__
        with open(newFolderName + "option.json", 'w') as f:
            json.dump(settingDic, f)


    if pathExists and canLoad:
        abunSplitter.mainFlow(newFolderName, newMummerLink)


    else:
        print "Sorry. The above folders or files are missing or options are not correct. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"

    print  "Time", time.time() - t0
    
    if houseKeeper.globalRunMPI == True:
        for i in range(1, nproc):
            data = "endall"
            comm.send(data, dest=i)


else:
    while True:
        data = comm.recv(source=0)

        if data == "endall":
            break

        elif len(data) > 0 and data[0] == "nucmerjob":
            mummerLink, folderName, outputName, referenceName, queryName, specialForRaw , specialName , refinedVersion= data[1:]
            alignerRobot.useMummerAlign(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw , specialName , refinedVersion)
            comm.send(data, dest=0)
        elif len(data) > 0 and data[0] == "gapjob":
            eachmatchpair,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename = data[1:]
            newdata = abunGraphLib.singleGapLookUp(eachmatchpair,folderName, N1,  mummerLink,  contigReadGraph, contigFilename,readsetFilename)
            comm.send(newdata, dest=0)
            

    