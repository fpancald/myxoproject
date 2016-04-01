#!/bin/csh
#$ -t 1-100:1
#$ -r y
#$ -q long
#$ -pe smp 15
#$ -N dss010


module load matlab/8.4
cd ..
matlab -nodisplay -nosplash -nojvm -r "diffsavesplit(${SGE_TASK_ID},0,1,0);exit"
##see diffsavesplit for description of input variables. order ntask, state, type, Rswitch.
### 010,0 fixed, half, noreceptor
