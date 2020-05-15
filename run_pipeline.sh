#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q workq -o /gpfs/home/natem/logs/snakelog.txt -j oe

cd /gpfs/home/natem/scripts/zika-pipeline/
module load gcc/6.3.0
snakemake -k -j 50 --cluster-config cluster.json --cluster "qsub -V -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}" --configfile snakemake_config.json
