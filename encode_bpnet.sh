#!/bin/bash -l

#$ -P siggers       # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00   # Specify the hard time limit for the job
#$ -l gpus=1       # Specify how many GPUs
#$ -l gpu_c=6.0       # Specify compute capacity of GPUs
#$ -pe omp 4       # Specify CPUs
#$ -N encode_h3k27ac           # Give job a name
#$ -j y               # Merge the error and output streams into a single file

cd /projectnb/siggers/thadz/siggers-rotation/

module load miniconda
module load samtools
module load bedtools
module load ucscutils

conda activate bpnet

samtools sort data/ENCFF600THN.bam -o data/ENCFF600THN_sorted.bam
bedtools genomecov -5 -bg -ibam data/ENCFF600THN_sorted.bam | LC_COLLATE=C sort -k1,1 -k2,2n > data/ENCFF600THN.bedGraph
bedGraphToBigWig data/ENCFF600THN.bedGraph data/GRCh38_EBV.chrom.sizes.tsv data/ENCFF600THN.bw

bpnet negatives -i data/ENCFF269RAB.bed -f /projectnb/siggers/genomes/hg38/sequences/genome.fa -b data/ENCFF600THN.bw -o data/ENCFF269RAB_matched_loci.bed -l 0.02 -w 2114 -v
bpnet fit -p jsons/bpnet_fit_encode.json
bpnet predict -p jsons/bpnet_predict_encode.json