#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler

PROJECT="infbench"

module load matlab/2017a
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)

export MATLABPATH=${MATLABPATH}:/${HOME}/${PROJECT}/matlab:${HOME}/MATLAB
source ${HOME}/MATLAB/setpath.sh

PROBLEMDIR="${HOME}/neurobench-problems"

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
        IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi
if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
        IID=${SLURM_ARRAY_TASK_ID}
fi

# Run the program

#PARAMS=$(awk "NR==${IID} {print;exit}" ${INPUTFILE})

cat<<EOF | matlab -nodisplay
addpath(genpath('${HOME}/vbmc/matlab'));
addpath(genpath('${HOME}/infbench/matlab'));
cd('${WORKDIR}');
id=$IID;
ii=mod(id-1,1000)+1;
Ns=$SAMPLES;
n=floor((id-1)/1000)+1;
infbench_akrami2018([],[],[ii,n,Ns]);
EOF
