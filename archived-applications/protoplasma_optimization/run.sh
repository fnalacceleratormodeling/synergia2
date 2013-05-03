# !/bin/bash

#q=amd32
#ppn=32
q=intel12
ppn=8
delta=0.0
offset=0.0
dpp=0.005
delta=1.0
#for proc in 32 #64 32 16 8 4 2 1
#do
#    synergia plasma.py createjob=1 submit=1 queue=$q \
#            jobdir=/data/cspark/results/nlopt/run \
#            numproc=$proc procspernode=$ppn
#done

proc=8
synergia fodo.py createjob=1 submit=1 queue=$q numproc=$proc \
         procspernode=$ppn jobdir=/data/cspark/results/nlopt/fodo
