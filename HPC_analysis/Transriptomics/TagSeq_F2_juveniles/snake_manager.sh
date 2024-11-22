#!/bin/bash
#SBATCH -t 250:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --mail-user=samuel.gurr@noaa.gov
#SBATCH --output=./"%x_out.%j"
#SBATCH --error=./"%x_err.%j"

source ~/.bashrc
mamba activate snakemake-7.32.4
snakemake --unlock # in case the job was canceled prior,
snakemake --cluster "sbatch -t {cluster.time} -N {cluster.nodes} --ntasks-per-node {cluster.ntasks-per-node} --mem {cluster.mem}" --latency-wait 20 --cluster-config config.yml -j 10
