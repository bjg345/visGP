#!/bin/bash

#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -pe make 1

module load conda_R
Rscript gen_finger_points.R

