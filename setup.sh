#!/bin/bash

###
git clone https://github.com/XingxingJian/script.git ~/snakemake
cd ~/snakemake
conda env create -f snakemake-tutorial.yaml
conda activate snakemake-tutorial

####download database
mkdir Database
cd Database
####hg38
mkdir hg38
mkdir k2_pluspf_20240605
cd hg38
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
####kraken2_db
cd ~/snakemake/Database/k2_pluspf_20240605
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz
gzip -d k2_pluspf_20240605.tar.gz

####
cd ~/snakemake

