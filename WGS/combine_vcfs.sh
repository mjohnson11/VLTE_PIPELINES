#!/bin/bash
#SBATCH -J combine_vcfs  #job name for array
#SBATCH -n 4                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 2-40:00              # Runtime in D-HH:MM
#SBATCH -p desai       # Partition to submit to
#SBATCH --mem=20000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/combine_vcfs_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/combine_vcfs_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

CHROMO=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ../accessory_files/chromo_intervals.txt)

OUTD="../../Output/WGS/"

# Getting sample nammodule load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/gatk/gatk-4.1.3.0

java -Xmx19g -XX:ParallelGCThreads=3 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar GenomicsDBImport -R ${OUTD}reference/w303_vlte.fasta --sample-name-map ../accessory_files/vcf_map.txt --genomicsdb-workspace-path ${OUTD}chromo_dbs/${CHROMO}_db -L ${CHROMO}

java -Xmx19g -XX:ParallelGCThreads=3 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar GenotypeGVCFs -R ${OUTD}reference/w303_vlte.fasta -V gendb://${OUTD}chromo_dbs/${CHROMO}_db -O ${OUTD}chromo_vcfs/${CHROMO}_final.vcf --heterozygosity 0.005 2> ${OUTD}logs/${CHROMO}_GenotypeGVCFs.log
