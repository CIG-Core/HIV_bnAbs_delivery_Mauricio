#!/bin/bash
#SBATCH -t 05-00:00
#SBATCH --mail-user=idm17@pitt.edu
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1
#SBATCH --mem=150g
#SBATCH --cpus-per-task=16
#SBATCH --job-name=s3Vs2Day28_cv

module load gcc/10.2.0
module load r/4.3.0

Rscript SLIDEcv.R /ix/cigcore/proj/mauricio/code/SLIDE_CV_s3Vs2Day28.yaml 