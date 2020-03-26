#!/bin/bash
module purge
module load matlab/2017a
BASEPATH="${HOME}/MATLAB"
source ${BASEPATH}/setpath.sh
export MATLABPATH=${MATLABPATH}:${HOME}/MATLAB
matlab -nodisplay
