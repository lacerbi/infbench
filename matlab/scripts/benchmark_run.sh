#!/bin/bash
PROJECT="infbench"
SHORTNAME=IB
BASEDIR="${HOME}/${PROJECT}/matlab"
SOURCEDIR="${BASEDIR}/"
JOBSCRIPT="${BASEDIR}/scripts/myjob.sh"

#Job parameters
RUN=${1}
INPUTFILE="${SCRATCH}/${PROJECT}/joblist-${1}.txt"
MAXID=$(sed -n $= ${INPUTFILE})

RUNTIME=6:00:00
MAXRT=NaN
VERBOSE=0
MAXFUNMULT="[]"

NODES="1"
PPN="1"
MEM="4000MB"
RESOURCES="nodes=${NODES}:ppn=${PPN},mem=${MEM},walltime=${RUNTIME}"
MAXSAMPLES="5000"
#MAXSAMPLES="100000"

#if [[ -z ${1} ]]; then
        JOBLIST="1-$MAXID"
        NEWJOB=1
#else
#        JOB=${1}
#        NEWJOB=0
#        echo "JOB=${JOB}" >> ${BASEDIR}/reruns.log
#fi

#Convert from commas to spaces
JOBLIST=${JOBLIST//,/ }
echo JOBS $JOBLIST

WORKDIR="${SCRATCH}/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}

JOBNAME=${SHORTNAME}${RUN}

#Submit jobs based on cluster name (defined in .bashrc)
if [ ${CLUSTER} = "Prince" ]; then
        sbatch --error=slurm-%A_%a.err --verbose --array=${JOBLIST} --mail-type=FAIL --mail-user=${USEREMAIL} --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},MATLAB_MODULE=${MATLAB_MODULE},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,MAXRT=$MAXRT,MAXSAMPLES=${MAXSAMPLES},INPUTFILE=${INPUTFILE},VERBOSE=${VERBOSE},MAXFUNMULT=${MAXFUNMULT} --job-name=${JOBNAME} ${JOBSCRIPT}
elif [ ${CLUSTER} = "Turso" ]; then
        sbatch --error=slurm-%A_%a.err --verbose --array=${JOBLIST} --clusters=ukko,vorna --partition=short,medium --mail-type=FAIL --mail-user=${USEREMAIL} --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},MATLAB_MODULE=${MATLAB_MODULE},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,MAXRT=$MAXRT,MAXSAMPLES=${MAXSAMPLES},INPUTFILE=${INPUTFILE},VERBOSE=${VERBOSE},MAXFUNMULT=${MAXFUNMULT} --job-name=${JOBNAME} ${JOBSCRIPT}
else
	qsub -t ${JOBLIST} -q normal -v PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,MAXRT=$MAXRT,MAXSAMPLES=${MAXSAMPLES},INPUTFILE=${INPUTFILE},VERBOSE=${VERBOSE},USEPRIOR=${USEPRIOR},TOLFUN=${TOLFUN},MAXFUNMULT=${MAXFUNMULT},STOPSUCCRUNS=${STOPSUCCRUNS} -l ${RESOURCES} -M ${USER}@nyu.edu -N ${JOBNAME} ${JOBSCRIPT}
fi
