#!/bin/bash

#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y

module load conda_R
Rscript write_grid.R

