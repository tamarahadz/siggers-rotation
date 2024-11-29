# Siggers Rotation Fall 2024: Training BPNet on Histone Mark Data
This repository walks through how to install BPNet onto the SCC, preprocess the data, train BPNet, and use the trained model to predict ChIP signal. The beginning of this tutorial assumes you have a configured your ```.condarc``` file (Step 1 in Emily's [Machine Learning in the SCC](https://github.com/ehk-kim/BUrotations/blob/main/Siggers/tutorials.md) tutorial). 

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


