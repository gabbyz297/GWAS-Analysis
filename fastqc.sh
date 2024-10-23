#!/bin/bash
#SBATCH -N 1
#SBATCH -n 17
#SBATCH -t 00:40:00
#SBATCH -p single
#SBATCH -A loni_chiguiro24
#SBATCH -o fastqc_pool1_102224.out
#SBATCH -e fastqc_pool1_102224.err

#This script will be run in parallel
module load parallel/20190222/intel-19.0.5

#Set the maximum heap size to 192 gigabytes
java -Xmx192g

#Specify where conda is installed
source /home/gabby297/miniconda3/etc/profile.d/conda.sh

#Activate conda environment
conda activate fastqc

#Specify where FastQC is installed; FastQC is installed locally you can skip all this and put the full path into the code on line 27 at the start of the quotes
/home/gabby297/miniconda3/envs/fastqc/bin/fastqc

#Change directory to where files are located
cd /work/gabby297/delaware_pool1/

#Use cat to call the list of samples parallel will use to call files and parallel --jobs 16 to run 16 jobs in parallel. Use {}.fastq.gz to specify the file type at the end of all the files and use -o to specify where the output files should be written.

cat sample_list.txt | parallel --jobs 16 "fastqc {}.fastq.gz -o /work/gabby297/delaware_pool1/fastqc/"
