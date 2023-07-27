#!/bin/bash

#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -pe local 2
#$ -j y
#$ -R y
#$ -t 1:1500
export _JAVA_OPTIONS=-Xmx3g

matlab -nodesktop -nosplash -nodisplay -r "run glgp.m"

