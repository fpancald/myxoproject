#!/bin/csh
#$ -t 1-100:1
#$ -r y
#$ -q long
#$ -pe smp 17
#$ -N diffsaverec


module load matlab
cd ..
matlab -nodisplay -nosplash -nojvm -r "diffusionsave2(${SGE_TASK_ID});exit"
