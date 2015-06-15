#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m ae
#$ -r y
#$ -q long
#$ -N job1myxo


module load matlab

matlab -nodisplay -nosplash < job1.m
