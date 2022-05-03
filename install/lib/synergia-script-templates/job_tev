#!/bin/sh
#PBS -o synergia.out
#PBS -e synergia.err
#PBS -l nodes=@@numnode@@
#PBS -l walltime=24:00:00
#PBS -q intel12
#PBS -A booster 

# This script is an appropriate template for running PBS jobs on the
# tev cluster.

echo PBS_NODEFILE is "$PBS_NODEFILE"
if [ -n "$PBS_NODEFILE" ]
then
    echo "begin nodes"
    cat $PBS_NODEFILE
    echo "end nodes"
fi
MPI_HOME=/usr/local/openmpi

export PATH=${MPI_HOME}/bin:$PATH
export PATH=/usr/local/python-2.7.0/bin:$PATH
export LD_LIBRARY_PATH=${MPI_HOME}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/python-2.7.0/lib:$LD_LIBRARY_PATH


. "/home/macridin/synergia2-crb4/setup.sh"

cd $PBS_O_WORKDIR
echo "job start at `date`"
mpirun --loadbalance -np @@numproc@@ @@synergia_executable@@ @@script@@ @@args@@ 
retval=$?
if [ x"$retval" != x"0" ]; then
	echo "job failed with return value $retval"
fi
echo "job end at `date`."



