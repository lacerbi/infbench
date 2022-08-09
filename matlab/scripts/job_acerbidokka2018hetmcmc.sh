#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

PROJECT="infbench"

# Use Intel compiler
if [ ${CLUSTER} = "Mercer" ]; then
        module load matlab gcc
        export LD_PRELOAD=$GCC_LIB/libstdc++.so
elif [ ${CLUSTER} = "Prince" ]; then
        module load matlab/2019a
        export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)
else
        module load ${MATLAB_MODULE}
        #export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)
fi
export MATLABPATH=${MATLABPATH}:/${HOME}/${PROJECT}/matlab:${HOME}/MATLAB:${WRKDIR}/MATLAB
#source ${HOME}/MATLAB/setpath.sh
export MATLABPATH="${MATLABPATH}:${HOME}/infbench/matlab:${HOME}/vbmc:${HOME}/vbmc-dev/acq:${HOME}/vbmc-dev/ent:${HOME}/vbmc-dev/misc:${HOME}/vbmc-dev/warp:${HOME}/infbench-dev/problems"

PROBLEMFOLDER="[]"

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
        IID=${PBS_ARRAYID}
elif [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
elif [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
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
infbench_acerbidokka2018het([],[],[ii,Ns]);
EOF
