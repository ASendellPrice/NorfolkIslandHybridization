#!/bin/bash

#If directory "data/assembly" does not exist then create it
if [[ ! -d data/assembly ]]
then
    mkdir data/assembly
fi

#Move into that directory
cd data/assembly

#Download the silvereye and zebra finch genome assemblies
#Note: we will be using the VGP zebra finch assembly which is assembled to chromosome-level and includes both Z and W
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.4_bTaeGut1.4.pri/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/735/GCA_001281735.1_ASM128173v1/GCA_001281735.1_ASM128173v1_genomic.fna.gz

#Set path to chromosemble tool and satsuma bin directory
CHROMOSEMBLE=/data/Users/BIN/satsuma2/bin/Chromosemble
export SATSUMA2_PATH=/data/Users/BIN/satsuma2/bin/

#Set query (silvereye) and target (zebra finch) assemblies
QUERY=GCA_001281735.1_ASM128173v1_genomic.fna.gz
TARGET=GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz

#Define output directory name
OUT=Zlat_2_Tgutt_ZW

#Run chromosemble
$CHROMOSEMBLE -t $TARGET -q $QUERY -o $OUT

