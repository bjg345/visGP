#!/bin/bash

#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y

export _JAVA_OPTIONS=-Xmx3g
module load conda_R
Rscript chesapeake.R

