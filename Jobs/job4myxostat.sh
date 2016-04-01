#!/bin/csh
#$ -M fpancald@nd.edu
#$ -m ae
#$ -r y
#$ -q long
#$ -pe smp 12
#$ -N job4


module load matlab

matlab -nodisplay -nosplash < job4.m
