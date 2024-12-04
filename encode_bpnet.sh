#!/bin/bash -l

#$ -P siggers       # Specify the SCC project name you want to use
#$ -m eas   # Send email when job finishes
#$ -o logs/encode_h3k27ac.log # Name for log file
#$ -l h_rt=12:00:00   # Specify the hard time limit for the job
#$ -l gpus=1       # Specify how many GPUs
#$ -l gpu_c=6.0       # Specify compute capacity of GPUs
#$ -pe omp 4       # Specify CPUs
#$ -N encode_h3k27ac           # Give job a name
#$ -j y               # Merge the error and output streams into a single file

cd /projectnb/siggers/thadz/siggers-rotation/
source ./common.sh
prerun

module load miniconda
module load samtools
module load bedtools
module load ucscutils

conda activate bpnet

run "samtools sort data/ENCFF600THN.bam -o data/ENCFF600THN_sorted.bam"
run "bedtools genomecov -5 -bg -ibam data/ENCFF600THN_sorted.bam | LC_COLLATE=C sort -k1,1 -k2,2n > data/ENCFF600THN.bedGraph"
run "bedGraphToBigWig data/ENCFF600THN.bedGraph data/GRCh38_EBV.chrom.sizes.tsv data/ENCFF600THN.bw"

run "bpnet negatives -i data/ENCFF269RAB.bed -f /projectnb/siggers/genomes/hg38/sequences/genome.fa -b data/ENCFF600THN.bw -o data/ENCFF269RAB_matched_loci.bed -l 0.02 -w 2114 -v"
run "bpnet fit -p jsons/bpnet_fit_encode.json"
run "predict -p jsons/bpnet_predict_encode.json"

run "python plot_logs.py -i h3k27ac_encode.log -o h3k27ac_encode.png"

postrun