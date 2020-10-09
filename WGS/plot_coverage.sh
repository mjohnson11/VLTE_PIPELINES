#!/bin/bash
#SBATCH -J plot_cov  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-05:00              # Runtime in D-HH:MM
#SBATCH -p desai       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/plot_cov_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/plot_cov_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent


source activate milo_py37

python plot_coverage.py ${SLURM_ARRAY_TASK_ID}
