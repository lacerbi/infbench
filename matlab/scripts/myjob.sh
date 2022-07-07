#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

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
#:${HOME}/vbmc/matlab/gplite:${HOME}/vbmc/matlab/utils"

#PROBLEMFOLDER="'${HOME}/neurobench-problems'"
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
PARAMS=$(awk "NR==${IID} {print;exit}" ${INPUTFILE})

echo ${PARAMS} ${VERBOSE} ${USEPRIOR}

cat<<EOF | matlab -nodisplay
%addpath(genpath('${HOME}/MATLAB'));
cd('${WORKDIR}');
options=struct('RootDirectory','${WORKDIR}','Display',${VERBOSE},'MaxFunEvalMultiplier',${MAXFUNMULT},'ProblemDirectory',${PROBLEMFOLDER},'MaxStoredSamples',${MAXSAMPLES});
${PARAMS}
infbench_run(${PARAMS},options);
EOF
