#!/bin/bash

#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1:1
module load conda_R
Rscript bora_csfeed.R $SGE_TASK_ID
