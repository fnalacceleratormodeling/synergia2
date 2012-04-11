# !/bin/bash

ppn=8
if [ $ppn = 8 ]; then
    q=intel12
    echo $q
else
    q=amd32
    echo $q
fi

maxturn=-1
numprocs=32
output_path=/data/cspark/results/mu2e
for solver in nosc #2d 3d
do
    dir=$output_path/checkpoint
    resume=$output_path/checkpoint.09
    synergia resumer.py createjob=1 submit=1 jobdir=$dir \
            verbosity=-1 \
            queue=$q \
            numproc=$numprocs \
            procspernode=$ppn \
            resume_dir=$resume \
            max_turns=$maxturn \
            walltime=00:30:00
done

