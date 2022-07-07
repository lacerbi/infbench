#!/bin/bash
PROJECT="infbench"
SHORTNAME=IB
BASEDIR="${HOME}/${PROJECT}/matlab"
SOURCEDIR="${BASEDIR}/"
JOBSCRIPT="${BASEDIR}/scripts/job_yachtmcmc.sh"

#Job parameters
RUNTIME=4:00:00
MAXRT=NaN
VERBOSE=0
WSAMPLES=200

NODES="1"
PPN="1"
MEM="4000MB"
RESOURCES="nodes=${NODES}:ppn=${PPN},mem=${MEM},walltime=${RUNTIME}"

RUN="yachtmcmc"

#if [[ -z ${1} ]]; then
        JOBLIST="$1"
        NEWJOB=1
#else
#        JOB=${1}
#        NEWJOB=0
#        echo "JOB=${JOB}" >> ${BASEDIR}/reruns.log
#fi

#Convert from commas to spaces
JOBLIST=${JOBLIST//,/ }
echo JOBS $JOBLIST

WORKDIR="${SCRATCH}/${PROJECT}/${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}

JOBNAME=${SHORTNAME}${RUN}

# running on Prince
sbatch --error=slurm-%A_%a.err --verbose --array=${JOBLIST} --mail-type=FAIL --mail-user=${USER}@nyu.edu --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},RUN=${RUN},WORKDIR=$WORKDIR,USER=$USER,WSAMPLES=$WSAMPLES --job-name=${JOBNAME} ${JOBSCRIPT}



