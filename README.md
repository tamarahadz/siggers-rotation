# Siggers Rotation Fall 2024: Training BPNet on H3K27ac ChIP-seq Data
This repository walks through how to install BPNet Lite onto the SCC, preprocess the data, train BPNet, and use the trained model to predict ChIP signal. The beginning of this tutorial assumes you have a configured your ```.condarc``` file (Step 1 in Emily's [Machine Learning in the SCC](https://github.com/ehk-kim/BUrotations/blob/main/Siggers/tutorials.md) tutorial). 

If you want to use this tool with the SCC On Demand interface, I highly reccommend reading through her full tutorial first. However, the rest of this tutorial assumes that you have connected to the SCC from your local computer with an [SSH client](https://www.bu.edu/tech/support/research/system-usage/connect-scc/ssh/). If you are unfamiliar with the SCC, I also recommend looking through [this workshop](https://programming-workshops.readthedocs.io/en/latest/workshops/01_linux_bash_scc/the_scc.html)

## Creating Your BPNet Conda Environment and Installing BPNet Lite
First, make sure you are in your project directory 
```console
bash-4.4~$ cd /projectnb/siggers/thadz
bash-4.4~$ pwd
/projectnb/siggers/thadz
```

Load the miniconda module 
```console
bash-4.4~$ module load miniconda
```

Create the conda environment using mamba. The bpnet-lite package requires Python >=3.9 and <3.13. Here I used 3.10. 
```console
bash-4.4~$ mamba create -y -n bpnet python=3.10
```

Once you have created your conda environment, you can activate it and install bpnet-lite into it using pip.
```console 
bash-4.4~$ conda activate bpnet
bash-4.4~$ pip install bpnet-lite
```

## Getting Your BPNet Input Files
The BPNet pipeline requires three input files - a FASTA file of the genome (the hg38 genome has already been downloaded under ```/projectnb/siggers/genomes/hg38/sequences/genome.fa```), BAM file of single-end or paired-end ChIP-seq reads which have already been aligned to the genome of interest, and a BED file of peaks (from the same experiment as the BAM file). 

Create a directory where you will store your data 
```console 
bash-4.4~$ mkdir data
```

Copy your BAM and BED files to that directory. For the ENCODE data, get files from their corresponding link on the ENCODE website 
```console 
bash-4.4~$ cd data
bash-4.4~$ wget https://www.encodeproject.org/files/ENCFF226FKB/@@download/ENCFF226FKB.bam
bash-4.4~$ wget https://www.encodeproject.org/files/ENCFF269RAB/@@download/ENCFF269RAB.bed.gz
bash-4.4~$ gunzip ENCFF269RAB.bed.gz
```
For the data generated in the Siggers lab, you will have BAM files and NarrowPeak files. The NarrowPeak files are in BED file format and can still be used in the BPNet pipeline. ChIP-seq experiment output can be found at ```/projectnb/siggers/Inge_ChIP/JK_ChIP_723/COFpendium/output/```. "UT" stands for "untreated" and "T" stands for "treated", and replicates are biological replicates.
```console
bash-4.4~$ cd data
bash-4.4~$ cp /projectnb/siggers/Inge_ChIP/JK_ChIP_723/COFpendium/output/peaks/H3K27ac_H3K27ac_UT1.narrowPeak .
bash-4.4~$ cp /projectnb/siggers/Inge_ChIP/JK_ChIP_723/COFpendium/output/bam/H3K27ac_H3K27ac_UT1_H3K27acU1.bam .
```

The actual pipeline requires us to convert the BAM file to a bigWig file of genome coverage before starting. To do this, we need one more file that specifies the sizes of every chromosome mentioned in the BAM file. 

For data from ENCODE, you can get this file from their database 
```console 
bash-4.4~$ cd data
bash-4.4~$ wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv
```

For data generated in the Siggers lab, we do not have this file on the SCC. However, you can download it from UCSC 
```console 
bash-4.4~$ cd data
bash-4.4~$ wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```
The Siggers lab uses a slightly different notation than UCSC to refer to chromosomes other than X, Y, and 1-21. For example, "chr15_KI270905v1_alt" is referred to as "KI270905.1" in the BAM files. To convert these names, you can use the following while loop in bash. For each chromosome, it replaces "v" with ".", removes "_random" and "_alt", then takes the second element of any chromosome name separated by an underscore, and writes the modified file to a new file named ```hg38_siggers.chrom.sizes```. 
```console 
bash-4.4~$ while IFS="" read -r p; do 
    echo "$p" | sed -e 's/v/./g' -e 's/_random//g' -e 's/_alt//g' | cut -d "_" -f 2 >> hg38_siggers.chrom.sizes; 
done < hg38.chrom.sizes
```
Now we are ready to begin preprocessing