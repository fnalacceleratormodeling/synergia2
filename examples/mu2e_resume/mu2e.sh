# !/bin/bash

ppn=8
if [ $ppn = 8 ]; then
    q=intel12
    #q=test_intel12
    echo $q
else
    q=amd32
    echo $q
fi

voltage=0.032 #0.01066667
turn=200
maxturn=100
ramping=50
chkperiod=20
tgrid=32
lgrid=128
np=1000
solver=3d
nu_h=9.65
nu_v=9.78
wire=0.016
steps=240
numprocs=32
output_path=/data/cspark/results/mu2e
dir=$output_path/checkpoint
for solver in nosc
do
    synergia mu2e.py createjob=1 submit=1 \
            jobdir=$dir \
            verbosity=1 \
            queue=$q \
            procspernode=$ppn \
            numproc=$numprocs \
            lattice_load=1 \
            turn_particles=1 \
            map_order=4 \
            spacecharge=$solver \
            num_steps=$steps \
            num_turns=$turn \
            max_turns=$maxturn \
            rampturns=$ramping \
            checkpointperiod=$chkperiod \
            macro_particles=$np \
            gridx=$tgrid gridy=$tgrid gridz=$lgrid \
            tuneh=$nu_h \
            tunev=$nu_v \
            wire_x=$wire \
            rf_voltage=$voltage \
            walltime=00:30:00
done

