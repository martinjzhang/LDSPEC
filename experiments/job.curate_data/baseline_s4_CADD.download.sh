#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH --array=1
#SBATCH -p short
#SBATCH --mem=16000
#SBATCH -o /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE



DATA_PATH=/n/scratch3/users/j/jz286/CADD
echo $DATA_PATH
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz.tbi -P $DATA_PATH 
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz -P $DATA_PATH 

