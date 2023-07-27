
#!/bin/bash

#$ -l mem_free=15G
#$ -l h_vmem=15G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 3001:4500

module load conda_R
Rscript nearest.R $SGE_TASK_ID

