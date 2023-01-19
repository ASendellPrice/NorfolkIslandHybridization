#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   02_filter_reads_novogene.sh
# Author:   Ash Sendell-Price
# Date:     18/01/2023
# System:   vettel (warwick / pushes results to s3 bucket)
# Summary:  Conducts initial QC of raw sequencing reads plus performs basic read filtering using
#           fastp. Detects / removes adapter contect and duplicated reads and trims 10bp. Filtered
#           reads are pushed to an s3 bucket for archiving.
#################################################################################################

#Load conda environment
mamba activate NorfolkIslandHybridisation

#Make required directories
mkdir novogene_filtered_reads
mkdir novogene_fastp_reports

#Initiate log file
date > novogene_filtered_reads/read_filtering.log

#For each of the novogene sequenced samples in resources/novogene_samples.txt
#do the following:
for SAMPLE in $(cat resources/novogene_samples.txt)
do
   
    #Make sample subdirectory in novogene_filtered_reads
    mkdir novogene_filtered_reads/${SAMPLE}
    
    #Make a temporary directory and move into it
    mkdir temp
    cd temp

    #Copy samples raw reads from S3 bucket to temp directory
    ~/aws-cli/bin/aws s3 sync s3://norfolkhybrids/novogene_raw_data/$SAMPLE/ ./

    #For each sample read pair do the following ...
    for ReadPair in $(ls ${SAMPLE}_*_1.fq.gz | cut -f1,2,3,4 -d'_')
    do
        #Use Fastp to conduct automated filtering of fastq files
        #Note: based on initial test we will trim the first 10bp from start of each read
        fastp \
        -i ${ReadPair}_1.fq.gz \
        -o ../novogene_filtered_reads/${SAMPLE}/Filtered_${ReadPair}_1.fq.gz \
        -I ${ReadPair}_2.fq.gz \
        -O ../novogene_filtered_reads/${SAMPLE}/Filtered_${ReadPair}_2.fq.gz \
        --trim_front1 10 \
        --trim_front2 10 \
        --detect_adapter_for_pe \
        --dedup
        
        #Rename / move fastp reports to fastp_log directory
        mv fastp.html ../novogene_fastp_reports/${ReadPair}.html
        mv fastp.json ../novogene_fastp_reports/${ReadPair}.fastp.json

        #Send fastp reports to s3 bucket for storage
        ~/aws-cli/bin/aws s3 cp ../novogene_fastp_reports/${ReadPair}.html s3://norfolkhybrids/novogene_fastp_reports/
        ~/aws-cli/bin/aws s3 cp ../novogene_fastp_reports/${ReadPair}.json s3://norfolkhybrids/novogene_fastp_reports/

        #Send filtered reads to s3 bucket
        ~/aws-cli/bin/aws s3 cp ../novogene_filtered_reads/${SAMPLE}/Filtered_${ReadPair}_1.fq.gz s3://norfolkhybrids/novogene_filtered_reads/${SAMPLE}/
        ~/aws-cli/bin/aws s3 cp ../novogene_filtered_reads/${SAMPLE}/Filtered_${ReadPair}_2.fq.gz s3://norfolkhybrids/novogene_filtered_reads/${SAMPLE}/

    done

    #Move out of temp directory and delete it
    cd ../
    rm -r temp

    #Add entry to log file
    echo "read filtering of sample" $SAMPLE "complete!" >> novogene_filtered_reads/read_filtering.log

done

#Summarise fastp reports using multiqc
cd novogene_fastp_reports
multiqc .
mv multiqc_report.html novogene_multiqc_report.html
~/aws-cli/bin/aws s3 cp novogene_multiqc_report.html s3://norfolkhybrids/novogene_fastp_reports/