#!/bin/bash
module purge
module load matlab/2014a
BASEPATH="${ROOTPATH}/MATLAB"
source ${BASEPATH}/setpath.sh
export MATLABPATH=${MATLABPATH}:${HOME}/MATLAB
matlab -nodisplay
