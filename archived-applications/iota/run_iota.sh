#!/bin/bash
echo "use like that:  ./run_iota.sh working_directory executable_name resume_executable lattice_xml_file input_options_file 
   resume_input_options_file numproc procspernode other_files(optional)"
if [ -z "$8" ]
  then
    echo "No enough arguments supplied"
    echo "use like that ./run_iota.sh working_directory executable_name resume_executable lattice_xml_file input_options_file 
       resume_input_options_file numproc procspernode other_files(optional)"
    exit
fi


working_directory=$1
executable=$2
resume_executable=$3
lattice_file=$4
input_options_file=$5
resume_input_options_file=$6
numproc=$7
procspernode=$8

numnodes=$((numproc/procspernode)) 
echo "numnodes=" $numnodes

suf=0

working_directory=$1"."$suf

while [ -d "$working_directory" ]
       do
           suf=$((suf+1))
           working_directory=$1"."$suf   
           echo "loop working directory=" $working_directory
       done
echo "create directory"
mkdir $working_directory
echo "present working directory=" $PWD
echo "loop working directory=" $working_directory
echo "executable=" $executable
echo "resume_executable"=$resume_executable
cp $executable $working_directory
cp $resume_executable $working_directory
#cp $lattice_file $working_directory/booster_lattice.xml
cp $lattice_file $working_directory/.
cp $input_options_file $working_directory/input_options
cp $resume_input_options_file $working_directory/input_resume_options

for (( i=9; i<=$#; i++ ))
do   
   if [ -f "${!i}" ] 
      then
        echo "copy file  ${!i} to  $working_directory"
        cp ${!i} $working_directory/.
   else
      echo " file ${!i} is not found"
   fi
done



echo "creating job script"
scriptname=$(basename "$working_directory")
scriptname+="_job"




MPI_HOME=/usr/local/openmpi

cat <<EOF >$working_directory/$scriptname
#!/bin/bash
#PBS -o synergia.out
#PBS -e synergia.err
#PBS -l nodes=$numnodes:amd32
#PBS -l walltime=24:00:00
#PBS -q amd32
#PBS -A booster

# This script is an appropriate template for running PBS jobs on the
# tev cluster.

echo PBS_NODEFILE is "\$PBS_NODEFILE"
if [ -n "\$PBS_NODEFILE" ]
then
    echo "begin nodes"
    cat \$PBS_NODEFILE
    echo "end nodes"
fi



export PATH=${MPI_HOME}/bin:$PATH



. "/home/macridin/synergia2-booster/setup.sh"

cd $working_directory

echo "job start at \`date\` "



$MPI_HOME/bin/mpirun -np $numproc -map-by core --bind-to core $working_directory/$executable

retval=\$?
if [ x"\$retval" != x"0" ]; then
        echo "job failed with return value $retval"
fi
echo "job end at \`date\`."

EOF

echo "creating resume script"
rscriptname=$(basename "$working_directory")
rscriptname+="_resume"
echo $rscriptname

cat <<EOF >$working_directory/$rscriptname
#!/bin/sh
#PBS -o rsynergia.out
#PBS -e rsynergia.err
#PBS -l nodes=$numnodes:amd32
#PBS -l walltime=24:00:00
#PBS -q amd32
#PBS -A booster


# This script is an appropriate template for running PBS jobs on the
# tev cluster.

echo PBS_NODEFILE is "\$PBS_NODEFILE"
if [ -n "\$PBS_NODEFILE" ]
then
    echo "begin nodes"
    cat \$PBS_NODEFILE
    echo "end nodes"
fi
MPI_HOME=/usr/local/openmpi

export PATH=${MPI_HOME}/bin:$PATH



. "/home/macridin/synergia2-booster/setup.sh"

cd $working_directory
echo "saving existing output files"
idtime=\$(date +"%m%d%Y_%H-%M-%S")
if [ -f rsynergia.out ]; then 
  cp rsynergia.out rsynergia_\$idtime.out
  echo "synergia.out file saved"
fi
if [ -f rsynergia.err ]; then 
  cp rsynergia.err rsynergia_\$idtime.err
  echo "synergia.err file saved"
fi


echo "job start at \`date\` "

mpirun -np $numproc -map-by core --bind-to core $working_directory/$resume_executable

retval=\$?
if [ x"\$retval" != x"0" ]; then
        echo "job failed with return value $retval"
fi
echo "job end at \`date\`."


EOF
