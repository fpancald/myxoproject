#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m abe
#$ -r y
#$ -q long
#$ -pe smp 17
#$ -N plotscript


module load matlab
cd ..
matlab -nodisplay -nosplash < plotscript.m
