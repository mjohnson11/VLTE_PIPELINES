#!/bin/bash
#SBATCH -J all_smoove  #job name for array
#SBATCH -n 1                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-04:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=16000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../Output/shell_outs/all_smoove.out      # File to which STDOUT will be written
#SBATCH -e ../../../Output/shell_outs/all_smoove.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

MYBAMS=""
for ((i=1;i<540;i++));
do
    BAMPATH=$(sed -n ${i}'{p;q}' ../../accessory_files/bam_map_combined_option.txt | awk '{print $2}')
    MYBAMS+=" ../${BAMPATH}"
done

singularity exec ~/smoove_latest.sif smoove call --outdir ../../../Output/WGS/smoove_output/ --name all_samples --fasta ../../../Output/WGS/reference/w303_vlte.fasta -p 1 -x --genotype ${MYBAMS}
