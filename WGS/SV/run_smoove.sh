#!/bin/bash
#SBATCH -J run_smoove  #job name for array
#SBATCH -n 2                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-01:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../Output/shell_outs/run_smoove_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../../Output/shell_outs/run_smoove_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

WELL=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ../../accessory_files/Wells.txt)

singularity exec ~/smoove_latest.sif smoove call --outdir ../../../Output/WGS/smoove_output/ --name ${WELL} --fasta ../../../Output/WGS/reference/w303_vlte.fasta -p 1 --genotype ../../../Output/WGS/work/*${WELL}.final.bam
