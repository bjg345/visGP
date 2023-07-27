#!/bin/bash

#$ -l mem_free=24G
#$ -l h_vmem=24G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 3001
module load conda_R
Rscript cs_fit.R $SGE_TASK_ID
