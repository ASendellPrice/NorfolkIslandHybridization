#!/bin/bash

#Load conda environment
mamba activate NorfolkIslandHybridisation

#Make directory for analysis and move into it
mkdir satsuma_chromosemble
cd satsuma_chromosemble

#Download assemblies
wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.4_bTaeGut1.4.pri/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz
wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/281/735/GCA_001281735.1_ASM128173v1/GCA_001281735.1_ASM128173v1_genomic.fna.gz

#Unzip the fasta files and replace lower case nucleotides with uppercase
#i.e. remove soft clipping
zcat GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz \
| awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
> GCA_003957565.4_bTaeGut1.4.pri_genomic.fna
zcat GCA_001281735.1_ASM128173v1_genomic.fna.gz \
| awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
> GCA_001281735.1_ASM128173v1_genomic.fna

#Remove old fastas
rm *fna.gz

#Run chromosemble
#Note: for some reason the -t and -q flags are mixed up in this version of Satsuma2
export SATSUMA2_PATH=~/miniconda3/envs/NorfolkIslandHybridisation/bin
Chromosemble \
-q GCA_001281735.1_ASM128173v1_genomic.fna \
-t GCA_003957565.4_bTaeGut1.4.pri_genomic.fna \
-o Zlat_2_Tgutt_ZW -pseudochr

#Update chromosome names:
#Chromosemble outputs long ass chromosome names that cause issues
#with downstream analyses, so we will replace these with short names
#i.e. chr1 chr2 etc.
cp Zlat_2_Tgutt_ZW/pseudochromosomes.fasta Zlat_2_Tgutt_ZW.fasta
cat Zlat_2_Tgutt_ZW.fasta | grep "chromosome" | sed 's/>//g' > long.chrom.names
for CHROM in $(cat long.chrom.names)
do
    #Extract chromosome name / number from long name
    CHROM_SHORT=$(echo $CHROM | cut -d "," -f 1 | cut -d "_" -f 7)
    #Use sed to replace long chrom name with short chrom name 
    #This may take several minutes
    sed -i "s/$CHROM/chr$CHROM_SHORT/g" Zlat_2_Tgutt_ZW.fasta
done

#Shorten scaffold names by removing ",_whole_genome_shotgun_sequence"
sed -i "s/,_whole_genome_shotgun_sequence//g" Zlat_2_Tgutt_ZW.fasta
sed -i "s/,_complete_sequence//g" Zlat_2_Tgutt_ZW.fasta

#Remove un-needed file
rm long.chrom.names

#Recompress fasta files
gzip GCA_003957565.4_bTaeGut1.4.pri_genomic.fna
gzip GCA_001281735.1_ASM128173v1_genomic.fna
gzip Zlat_2_Tgutt_ZW.fasta
