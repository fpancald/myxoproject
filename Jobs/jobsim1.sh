#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m ae
#$ -r y
#$ -q long
#$ -pe smp 12
#$ -N jobsim1


module load matlab

matlab -nodisplay -nosplash < jobsim1.m
