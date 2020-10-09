#!/bin/bash
#SBATCH -J ann_vcfs_lumpy  #job name for array
#SBATCH -n 2                   # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-01:05              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=3000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../Output/shell_outs/annotate_vcfs_lumpy.out      # File to which STDOUT will be written
#SBATCH -e ../../../Output/shell_outs/annotate_vcfs_lumpy.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

# Getting sample nammodule load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01

OUTD="../../../Output/WGS/smoove_output/"

java -Xmx4g -jar ~/snpEff/snpEff.jar -c ~/snpEff/snpEff.config -v w303_vlte ${OUTD}all_samples.smoove.square.vcf > ${OUTD}all_samples.smoove.square.ann.vcf
java -Xmx4g -jar ~/gatk/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantsToTable -V ${OUTD}all_samples.smoove.square.ann.vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F SVTYPE -F SVLEN -F END -F ID -GF RO -GF AO -O ${OUTD}all_samples.smoove.square.ann.tsv


