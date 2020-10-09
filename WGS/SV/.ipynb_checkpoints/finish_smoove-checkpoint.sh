#!/bin/bash
#SBATCH -J finish_smoove  #job name for array
#SBATCH -n 2                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-01:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=16000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../Output/shell_outs/finish_smoove.out      # File to which STDOUT will be written
#SBATCH -e ../../../Output/shell_outs/finish_smoove.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent


singularity exec ~/smoove_latest.sif smoove paste --name all_samples --outdir ../../../Output/WGS/smoove_output/ ../../../Output/WGS/smoove_output/genotyped/*.vcf.gz

