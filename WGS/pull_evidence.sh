#!/bin/bash
#SBATCH -J pull_evidence  #job name for array
#SBATCH -n 2                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-05:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=16000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/pull_evidence_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/pull_evidence_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load samtools

# Pulls lines from the original SAM that might be evidence for a particular mutation in the filtered set (per well)
# this is lenient, not all of these reads will actually be used, but that's ok, just showing everything close is fine
python make_evidence_sam_scripts.py ${SLURM_ARRAY_TASK_ID}
