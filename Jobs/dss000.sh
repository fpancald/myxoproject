#!/bin/csh
#$ -t 1-100:1
#$ -r y
#$ -q long
#$ -pe smp 15
#$ -N dss000


module load matlab/8.4
cd ..
matlab -nodisplay -nosplash -nojvm -r "diffsavesplit(${SGE_TASK_ID},0,0,0);exit"
##see diffsavesplit for description of input variables. order ntask, state, type, Rswitch.
### 0,0,0 fixed, uni, noreceptor
