#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   02_filter_reads_novogene.sh
# Author:   Ash Sendell-Price
# Date:     18/01/2023
# System:   vettel (warwick)
# Summary:  ADD
#################################################################################################

#Load conda environment
mamba activate NorfolkIslandHybridisation

#Create directories to store sample bam files
mkdir mapping
mkdir mapping/Zlat_2_Tgutt_ZW

#Index pseudochromosome assembly using bwa-mem2
#bwa-mem2 index assemblies/Zlat_2_Tgutt_ZW.fasta.gz

#Initiate log file
date > mapping/Zlat_2_Tgutt_ZW/read_mapping.log

#For each of the novogene sequenced samples in resources/novogene_samples.txt
#do the following:
for SAMPLE in $(cat resources/novogene_samples.txt)
do
    #Make a temporary directory and move into it
    mkdir mapping_temp
    cd mapping_temp

    #For each read pair do the following ...
    for ReadPair in $(ls ../novogene_filtered_reads/${SAMPLE}/Filtered_${SAMPLE}_*_1.fq.gz | rev | cut -d "_" -f2- | rev)
    do
        BASE_NAME=$(basename $ReadPair)
        bwa-mem2 mem \
        ../assemblies/Zlat_2_Tgutt_ZW.fasta.gz \
        ${ReadPair}_1.fq.gz \
        ${ReadPair}_2.fq.gz \
        -R "@RG\tID:${BASE_NAME}\tSM:${SAMPLE}" \
        | samtools view -bS - \
        | samtools sort - > ${BASE_NAME}.bam
    done    

    #If number of BAMs is greater than 1 then merge bams into a single file else
    #just rename the bam file
    if [ $(ls Filtered_${SAMPLE}*.bam | wc -l) -gt 1 ]
    then
        samtools merge ${SAMPLE}.bam Filtered_${SAMPLE}*.bam -f
    else
        mv Filtered_${SAMPLE}*.bam ${SAMPLE}.bam
    fi

    #Sort bam file and remove unsorted version
    samtools sort ${SAMPLE}.bam > ../mapping/Zlat_2_Tgutt_ZW/${SAMPLE}.sorted.bam

    #Index bam file
    samtools index ../mapping/Zlat_2_Tgutt_ZW/${SAMPLE}.sorted.bam

    #Transfer full bam and index to s3 bucket
    ~/aws-cli/bin/aws s3 cp ../mapping/Zlat_2_Tgutt_ZW/${SAMPLE}.sorted.bam s3://norfolkhybrids/mapping/Zlat_2_Tgutt_ZW/
    ~/aws-cli/bin/aws s3 cp ../mapping/Zlat_2_Tgutt_ZW/${SAMPLE}.sorted.bam.bai s3://norfolkhybrids/mapping/Zlat_2_Tgutt_ZW/

    #Tidy up
    cd ../
    rm -r mapping_temp
    
    #Add entry to log file
    echo "mapping completed for sample" $SAMPLE "!" >> mapping/Zlat_2_Tgutt_ZW/read_mapping.log

done
