#!/bin/bash
#SBATCH -J make_combined_bam_and_pilon  #job name for array
#SBATCH -n 1     # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-15:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=67000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/make_combined_bam_and_pilon.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/make_combined_bam_and_pilon.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load samtools
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01
module load jdk/1.8.0_45-fasrc01

OUTD="../../Output/WGS/reference_refinement/" 

samtools merge ${OUTD}combined_ref.bam ${OUTD}work/*final.bam 

java -Xmx20g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam I=${OUTD}combined_ref.bam O=${OUTD}combined_ref.final.bam SORT_ORDER=coordinate CREATE_INDEX=true 2> ${OUTD}logs/combined_ref_sorting.log

java -Xmx20g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar ValidateSamFile I=${OUTD}combined_ref.final.bam O=${OUTD}combined_ref.validate.txt MODE=SUMMARY 2> ${OUTD}logs/combined_ref_validate.log

java -Xmx66g -jar ~/pilon/pilon-1.23.jar --genome ../accessory_files/orig_w303_ref/w303_ref.fasta --bam ${OUTD}combined_ref.final.bam --output w303_vlte --outdir ${OUTD}pilon_output/ --changes --tracks

python update_gff_and_fix_fasta_headers.py
