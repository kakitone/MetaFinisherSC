#!/usr/bin/env python
import numpy
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

comm = MPI.COMM_WORLD
me = comm.Get_rank()
nproc = comm.Get_size()

'''
@ runnit 
module load python
module load mpi4py
aprun -n 24 python-mpi ./mpiStarter.py

@ Commandlnie 
qsub runnit

@ Local environment
mpiexec  -n 5 python  srcRefactor/repeatPhaserLib/mpiStarter.py 
'''


randNum = numpy.zeros(1)

if me == 0 : 
	print "I am the master"
	# Block until other jobs are all done
	total = 0
	for i in range(1, nproc):
		comm.Recv(randNum, ANY_SOURCE)
		total += randNum[0]

	print "This is master : ", total

	### Run usual routines 
	

	


else: 
	
	### Run each of the alignment steps. 
	
	
	print "me, nproc", me, nproc
	# Perform alignment according to the , when done, 
	# send a signal to the master
	randNum = numpy.random.random_sample(1)
	print "Process", me, "drew the number", randNum[0]
	comm.Send(randNum, dest=0)

