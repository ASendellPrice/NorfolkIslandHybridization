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

#Unzip the fasta files
gunzip GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz
gunzip GCA_001281735.1_ASM128173v1_genomic.fna.gz

#Set path to chromosemble tool and satsuma bin directory
CHROMOSEMBLE=/data/Users/BIN/satsuma2/bin/Chromosemble
export SATSUMA2_PATH=/data/Users/BIN/satsuma2/bin/

#Set query (silvereye) and target (zebra finch) assemblies
QUERY=GCA_001281735.1_ASM128173v1_genomic.fna
TARGET=GCA_003957565.4_bTaeGut1.4.pri_genomic.fna

#Define output directory name
OUT=Zlat_2_Tgutt_ZW

#Run chromosemble
$CHROMOSEMBLE -t $TARGET -q $QUERY -o $OUT

#Update chromosome names 
cp ${OUT}/pseudochromosomes.fasta ${OUT}.fasta
cat ${OUT}.fasta | grep "chromosome" | sed 's/>//g' > long.chrom.names
for CHROM in $(cat long.chrom.names)
do
    #Extract chromosome name / number from long name
    CHROM_SHORT=$(echo $CHROM | cut -d "," -f 1 | cut -d "_" -f 7)
    #Use sed to replace long chrom name with short chrom name 
    #This may take several minutes
    sed -i "s/$CHROM/chr$CHROM_SHORT/g" ${OUT}.fasta
done

#Remove un-needed file
rm long.chrom.names

