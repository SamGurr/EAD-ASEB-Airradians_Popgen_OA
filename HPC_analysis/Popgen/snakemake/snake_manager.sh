#!/bin/bash
#SBATCH -t 250:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output=./"%x_out.%j"
#SBATCH --error=./"%x_err.%j"

source ~/.bashrc
mamba activate snakemake-7.32.4
snakemake --unlock # in case the job was canceled prior,
# Chris P. meeting on 9/20 using pip executor to call higher mem partition, can call my own venv snakemake 8.???? and should work 
# BUT requires threads and resrouces inaddtion to input and output for each rule in the Snakefile

snakemake --cluster "sbatch -t {cluster.time} -N {cluster.nodes} --ntasks-per-node {cluster.ntasks-per-node} --mem {cluster.mem}" --latency-wait 60 --cores 8 --cluster-config config.yml -j 10 
