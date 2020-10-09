#!/bin/bash
#SBATCH -J add_depth  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00              # Runtime in D-HH:MM
#SBATCH -p desai       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../Output/shell_outs/wgs_depth_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../Output/shell_outs/wgs_depth_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent


source activate milo_py37
OUTD="../../Output/WGS/coverage/"

for ((i=$((${SLURM_ARRAY_TASK_ID}*10+1));i<$((${SLURM_ARRAY_TASK_ID}*10+11));i++));
do
  SAMP=$(sed -n ${i}'{p;q}' ../accessory_files/bam_map_combined_option.txt | awk '{print $1}')
  BAMPATH=$(sed -n ${i}'{p;q}' ../accessory_files/bam_map_combined_option.txt | awk '{print $2}')
  echo $i $SAMP $BAMPATH
  samtools depth ${BAMPATH} > ${OUTD}${SAMP}_depth.tsv
done
