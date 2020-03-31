#!/bin/bash
#SBATCH -J vcf_parse  #job name for array
#SBATCH -n 1                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=16000               # Memory pool for all cores (see also --mem-per-cpu)
#aSBATCH -o vcf_parse.out      # File to which STDOUT will be written
#aSBATCH -e vcf_parse.err      # File to which STDERR will be written
#SBATCH -o ../../Output/shell_outs/vcf_parse.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/vcf_parse.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

#combines all the separate chromosome ones, adds some annotations, filters and splits into per well output
python final_vcf_parsing.py  
