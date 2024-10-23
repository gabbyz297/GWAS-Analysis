#!/bin/bash
#SBATCH -N 1
#SBATCH -n 17
#SBATCH -t 00:40:00
#SBATCH -p single
#SBATCH -A loni_chiguiro24
#SBATCH -o fastqc_pool1_102224.out
#SBATCH -e fastqc_pool1_102224.err


module load parallel/20190222/intel-19.0.5

java -Xmx192g

source /home/gabby297/miniconda3/etc/profile.d/conda.sh

conda activate fastqc

/home/gabby297/miniconda3/envs/fastqc/bin/fastqc

cd /work/gabby297/delaware_pool1/

cat sample_list.txt | parallel --jobs 16 "fastqc {}.fastq.gz -o /work/gabby297/delaware_pool1/fastqc/"
