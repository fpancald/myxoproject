#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m ae
#$ -r y
#$ -q long
#$ -pe smp 12
#$ -N jobdiff2


module load matlab

matlab -nodisplay -nosplash < jobdiff2.m
