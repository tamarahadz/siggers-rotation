#!/bin/bash -l

#$ -P siggers       # Specify the SCC project name you want to use
#$ -m eas   # Send email when job finishes
#$ -o logs/siggers_h3k27ac.log # Name for log file
#$ -l h_rt=12:00:00   # Specify the hard time limit for the job
#$ -l gpus=1       # Specify how many GPUs
#$ -l gpu_c=6.0       # Specify compute capacity of GPUs
#$ -pe omp 4       # Specify CPUs
#$ -N siggers_h3k27ac           # Give job a name
#$ -j y               # Merge the error and output streams into a single file

cd /projectnb/siggers/thadz/siggers-rotation/
source ./common.sh
prerun

module load miniconda
module load samtools
module load bedtools
module load ucscutils

conda activate bpnet

run "samtools sort data/H3K27ac_H3K27ac_UT1_H3K27acU1.bam -o data/H3K27ac_H3K27ac_UT1_H3K27acU1_sorted.bam"
run "bedtools genomecov -5 -bg -ibam data/H3K27ac_H3K27ac_UT1_H3K27acU1_sorted.bam | LC_COLLATE=C sort -k1,1 -k2,2n > data/H3K27ac_H3K27ac_UT1_H3K27acU1.bedGraph"
run "bedGraphToBigWig data/H3K27ac_H3K27ac_UT1_H3K27acU1.bedGraph data/hg38_siggers.chrom.sizes data/H3K27ac_H3K27ac_UT1_H3K27acU1.bw"

run "bpnet negatives -i data/H3K27ac_H3K27ac_UT1.narrowPeak -f /projectnb/siggers/genomes/hg38/sequences/genome.fa -b data/H3K27ac_H3K27ac_UT1_H3K27acU1.bw -o data/H3K27ac_H3K27ac_UT1_H3K27acU1_matched_loci.bed -l 0.02 -w 2114 -v"
run "bpnet fit -p jsons/bpnet_fit_siggers.json"
run "bpnet predict -p jsons/bpnet_predict_siggers.json"

run "python plot_logs.py -i h3k27ac_siggers.log -o h3k27ac_siggers.png"

postrun