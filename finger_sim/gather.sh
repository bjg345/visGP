#!/bin/bash

#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y

module load conda_R
Rscript gather.R
