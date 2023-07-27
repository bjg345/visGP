#!/bin/bash

#$ -l mem_free=40G
#$ -l h_vmem=40G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y

module load conda_R
Rscript cov_comp.R $SGE_TASK_ID

