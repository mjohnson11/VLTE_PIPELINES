#!/bin/bash
#SBATCH -J wgs_pipe_setup  #job name for array
#SBATCH -n 2                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-00:20              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=4000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../shell_outputs/wgs_pipe_setup.out      # File to which STDOUT will be written
#SBATCH -e ../shell_outputs/wgs_pipe_setup.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

OUTD="../../Output/WGS/"

module load samtools

samtools faidx ${OUTD}reference/w303_vlte.fasta

module load bwa/0.7.15-fasrc02

# indexing the ref. genome with BWA needs to be done once
# this will create 5 additional files

bwa index -p ${OUTD}reference/w303_vlte ${OUTD}reference/w303_vlte.fasta

module load jdk/1.8.0_45-fasrc01
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01

java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=${OUTD}reference/w303_vlte.fasta O=${OUTD}reference/w303_vlte.dict
