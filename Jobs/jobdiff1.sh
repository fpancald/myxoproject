#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m ae
#$ -r y
#$ -q long
#$ -pe smp 12
#$ -N jobdiff1


module load matlab

matlab -nodisplay -nosplash < jobdiff1.m
