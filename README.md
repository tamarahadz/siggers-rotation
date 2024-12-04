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

Copy your BAM and BED files to that directory. For the ENCODE data, get files from their corresponding link on the ENCODE website. For example, for the dataset [ENCSR000AKY](https://www.encodeproject.org/experiments/ENCSR000AKY/), the experiment has three biological replicates and has Control BAMs ([rep1](https://www.encodeproject.org/files/ENCFF226FKB/@@download/ENCFF226FKB.bam),[rep2](https://www.encodeproject.org/files/ENCFF666YXZ/@@download/ENCFF666YXZ.bam),[rep3](https://www.encodeproject.org/files/ENCFF420ARY/@@download/ENCFF420ARY.bam)), Experimental BAMs ([rep1](https://www.encodeproject.org/files/ENCFF600THN/@@download/ENCFF600THN.bam),[rep2](https://www.encodeproject.org/files/ENCFF704LGA/@@download/ENCFF704LGA.bam),[rep3](https://www.encodeproject.org/files/ENCFF232RQF/@@download/ENCFF232RQF.bam)), and BED files ([rep1](https://www.encodeproject.org/files/ENCFF269RAB/@@download/ENCFF269RAB.bed.gz),[rep2](https://www.encodeproject.org/files/ENCFF546CSD/@@download/ENCFF546CSD.bed.gz),[3](https://www.encodeproject.org/files/ENCFF652VDP/@@download/ENCFF652VDP.bed.gz)) 

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

## Preprocessing
After you have all of your files, you can create a bash script to do all of the following steps in this tutorial. You can follow the example in this repo, ```bpnet.sh```, which uses the data coming from the Siggers lab. You can submit your script to the SCC with 

```console 
bash-4.4~$ qsub bpnet.sh
```
You can also do this tutorial in interactive mode by following the rest of the steps in your console. 

Load the modules needed to convert the BAM file to a bigWig file from the SCC. 
```console 
bash-4.4~$ module load samtools
bash-4.4~$ module load bedtools
bash-4.4~$ module load ucscutils
```

Use samtools to sort the BAM file. This will create lots of temporary files first. If your BAM file is particularly large, you may not be able to run this command in interactive mode. 
```console 
bash-4.4~$ samtools sort bamfile.bam -o bamfile_sorted.bam
```

Use bedtools to convert the sorted BAM file to a bedGraph file. 

```console 
bash-4.4~$ bedtools genomecov -5 -bg -ibam bamfile_sorted.bam | LC_COLLATE=C sort -k1,1 -k2,2n > bamfile.bedGraph
```

Finally, use bedGraphToBigWig to convert the bedGraph to a bigWig file. This is where you need the chromosome sizes file (which you downloaded from the ENCODE database, or modified from the UCSC file)

```console 
bash-4.4~$ bedGraphToBigWig bamfile.bedGraph hg38_siggers.chrom.sizes bamfile.bw
```

The last step of preprocessing is the first step of the BPNet pipeline (extracting GC-matched negatives).
```console 
bash-4.4~$ bpnet negatives -i bedfile.bed -f genome.fa -b bamfile.bw -o matched_loci.bed -l 0.02 -w 2114 -v
```

Now we are ready to train BPNet.

## Training and Using BPNet
The bpnet command line tool takes arguments through a JSON file. Example JSONS for all steps in the pipeline can be found on the BPNet [GitHub Repo](https://github.com/jmschrei/bpnet-lite/tree/master/example_jsons). Be aware that some of the JSONS in this repository have typos (missing commas, etc.). You can create your json following these templates, or the ones in this repository. 

To use BPNet in interactive mode, you need to log into a compute node from the login node. These are the same parameters that are specified in the qsub script. You need to request a GPU with a compute capacity of at least 6. Along with it, we request 4 CPUs. These options are described on the BU SCC [TechWeb](https://www.bu.edu/tech/support/research/software-and-programming/gpu-computing/) website. 
```console 
bash-4.4~$ qrsh -l gpus=1 -l gpu_c=6.0 -pe omp 4
```

When you log into the compute node, it will automatically log you into your home directory and modules you loaded in the login node will not be loaded. Go back into your project directory and reactivate your conda environment 
```console 
bash-4.4~$ cd ../../../projectnb/siggers/thadz
bash-4.4~$ module load miniconda
bash-4.4~$ conda activate bpnet
```

Now fit the model using the parameters specified in your JSON files. 
```console 
bash-4.4~$ bpnet fit -p jsons/bpnet_fit_example.json
```

The torch file is the trained BPNet model, which you supply in the JSON for the predict step. I would encourage looking at the log file to examine the mean-squared error and Pearson correlation coefficient for testing before going to the prediction step. Run the prediction step with
```console 
bash-4.4~$ bpnet predict -p jsons/bpnet_predict_example.json
```
This will output two npz tiles, one that predicts the counts and one that prevents the profile of the ChIP.

## Future Steps
**Running with control files** - I was unable to figure out how to run the package with the control bigwig files. I used the same methodology used to process the experimental BAMs to process the control BAMs. Then in the JSON parameters, I specified the controls in the same way as I did with the experimentals. However, I would get the following type of pytorch error 

```console 
RuntimeError: Given groups=1, weight of size [1, 66, 75], expected input[64, 65, 2114] to have 66 channels, but got 65 channels instead
```

According to the original paper, the control track is supposed to be "input for ChIP–seq and PATCH–CAP for ChIP–nexus." I found this a bit unclear for our ChIP-seq data so I tried using a different experiment (ex: for the Siggers lab data, the EP300 experiment) or the files listed as ENCODE as controls (ex: ```ENCFF226FKB.bam``` is the control for ```ENCFF600THN.bam```), and both of these still gave me this error. 

**Interpreting results** - I have a vague idea about what the counts and profile refer to; however I was not able to go too into depth into understanding how these predictions can be interpreted biologically

**Using the full bersion of BPNet** - I was advised to use the lite version of BPNet, since the full version is built in TensorFlow and can often be hard to install. However, on the BPNet Lite repo, it states "This package is meant to be used for quickly exploring data sets using BPNet or ChromBPNet and as a springboard for prototyping new ideas that involve modifying the code." I have not yet looked into performance of BPNet vs BPNet Lite, but the main question I have is whether this means that BPNet Lite can't be used with DeepLIFT or TFModisco or to get biological insights about the data.

**Data performance** - The file system that hosts already preprocessed files mentioned in the paper can be found [here](https://mitra.stanford.edu/kundaje/avsec/chipnexus/paper/data/). There seem to be some differences between model performance on their ChIP-seq data and ChIP-seq data from ENCODE and ChIP-seq from the Siggers lab. ENCODE data and Siggers lab data both perform similarly to each other, but BPNet with and without controls has much higher Pearson correlation in predicting the profile. I am wondering whether this is due to data processing - on their repo it says that the alignment was done with ChroMAP - could this be affecting results?