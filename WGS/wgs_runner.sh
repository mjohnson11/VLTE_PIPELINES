#!/bin/bash
#SBATCH -J wgs_pipe  #job name for array
#SBATCH -n 2                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-10:00              # Runtime in D-HH:MM
#SBATCH -p desai       # Partition to submit to
#SBATCH --mem=8000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/wgs_pipe_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/wgs_pipe_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

# Getting sample name
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ../accessory_files/Samples.txt)
SUF=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ../accessory_files/Sample_Suffixes.txt)
PRE=$(sed -n ${SLURM_ARRAY_TASK_ID}'{p;q}' ../accessory_files/Sample_Prefixes.txt)

DATAB="../../Data/WGS/"
OUTD="../../Output/WGS/"

# 1 TRIMMING
# NGmerge (only good for paired-end data)

module load NGmerge/0.2-fasrc01

NGmerge -1 ${DATAB}raw_reads/${PRE}_${SAMP}${SUF}_R1_001.fastq.gz -2 ${DATAB}raw_reads/${PRE}_${SAMP}${SUF}_R2_001.fastq.gz -a -v -o ${OUTD}work/${SAMP}.trimmed

# 2 ALIGNMENT AGAINST A REFERENCE GENOME
# BWA: more suitable than Bowtie2 when GATK is used downstream  # because both BWA and GATK were created at Broad (same people)

module load bwa/0.7.15-fasrc02

# mem is the BWA function for alignment

bwa mem -M -t 1 -R "@RG\tID:HJ25MBGXC\tSM:${SAMP}\tPL:ILLUMINA" ${OUTD}reference/w303_vlte ${OUTD}work/${SAMP}.trimmed_1.fastq.gz ${OUTD}work/${SAMP}.trimmed_2.fastq.gz > ${OUTD}work/${SAMP}.sam 2> ${OUTD}logs/${SAMP}.bwa.log

# 3 CONVERSION FROM SAM TO BAM FILE
# can be done with Samtools or Picard
# here Picard is used and data are also sorted and indexed
# Picard runs on java

module load jdk/1.8.0_45-fasrc01
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01

java -Xmx7g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam I=${OUTD}work/${SAMP}.sam O=${OUTD}work/${SAMP}.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true

# 4 MARKING DUPLICATES
# also performed with Picard

java -Xmx7g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar MarkDuplicates I=${OUTD}work/${SAMP}.sorted.bam O=${OUTD}work/${SAMP}.dedup.bam METRICS_FILE=${OUTD}work/${SAMP}.dedup_metrics.txt REMOVE_DUPLICATES=false TAGGING_POLICY=All 2> ${OUTD}logs/${SAMP}_dedup.log

# 5 RESORTING AND REINDEXING
# SortSam (Picard) needs to be run again

java -Xmx7g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam I=${OUTD}work/${SAMP}.dedup.bam O=${OUTD}work/${SAMP}.final.bam SORT_ORDER=coordinate CREATE_INDEX=true 2> ${OUTD}logs/${SAMP}_final_sorting.log

# 6 VALIDATING THE BAM FILES
# also performed with Picard

java -Xmx7g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar ValidateSamFile I=${OUTD}work/${SAMP}.final.bam O=${OUTD}work/${SAMP}.validate.txt MODE=SUMMARY 2> ${OUTD}logs/${SAMP}_validate.log

# 7 IF "no errors found", ALL PREVIOUS SAM AND BAM FILES CAN NOW
# BE DELETED (IN PRINCIPLE IT IS ALWAYS POSSIBLE TO REGENERATE
# EVEN THE ORIGINAL FASTQ FILES FROM THE FINAL BAM FILE)

# 8 HAPLOTYPE CALLING
# we will use GATK4, the most recent version (Gold standard)
# 3 steps are needed: haplotype calling, database creation
# and genotyping (this is 1 step more than with GATK3)
# GATK4 needs to be copied to a known directory

module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
GATK4_HOME=~/gatk/gatk-4.1.3.0

# 8A generating the g.vcf files for the cohort samples

java -Xmx7g -XX:ParallelGCThreads=1 -jar $GATK4_HOME/gatk-package-4.1.3.0-local.jar HaplotypeCaller -R ${OUTD}reference/w303_vlte.fasta -I ${OUTD}work/${SAMP}.final.bam -O ${OUTD}work/${SAMP}.g.vcf --emit-ref-confidence GVCF
