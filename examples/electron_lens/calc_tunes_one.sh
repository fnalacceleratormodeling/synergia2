#PBS -o calc_tunes.out
#PBS -e calc_tunes.err
#PBS -l walltime=24:00:00
#PBS -A mi

if [ ${PBS_O_QUEUE} = "intel12" ] ; then
    ppn=12
elif [ ${PBS_O_QUEUE} = "amd32" ] ; then
    ppn=32
else
    echo "unknown queue name: ->${PBS_O_QUEUE}<-"
    exit
fi

nnodes=`wc -l ${PBS_NODEFILE} | cut -d' ' -f1`
nproc=$(($nnodes*$ppn))

echo "running on $nproc processors"

. ~/syn2-dev/setup.sh

cd ${PBS_O_WORKDIR}

mpirun -np ${nproc} synergia ~/electron_lens/calc_tunes_onefile.py bulk_track_0.h5
