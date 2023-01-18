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
mkdir filtered_reads
mkdir fastp_logs

#Initiate log file
date > filtered_reads/read_filtering.log

#For each of the novogene sequenced samples in resources/novogene_samples.txt
#do the following:

for SAMPLE in $(cat resources/novogene_samples.txt)
do
   
    #Make sample subdirectory in filtered_reads
    mkdir filtered_reads/${SAMPLE}
    
    #Make a temporary directory and move into it
    mkdir temp
    cd temp

    #Copy samples raw reads from S3 bucket to temp directory
    ~/aws-cli/bin/aws s3 sync s3://norfolkhybrids/raw_data/$SAMPLE/ ./

    #For each sample read pair do the following ...
    for ReadPair in $(ls ${SAMPLE}_*_1.fq.gz | cut -f1,2,3,4 -d'_')
    do
        #Use Fastp to conduct automated filtering of fastq files
        #Note: based on initial test we will trim the first 10bp from start of each read
        fastp \
        -i ${ReadPair}_1.fq.gz \
        -o ../filtered_reads/${SAMPLE}/Filtered_${ReadPair}_1.fq.gz \
        -I ${ReadPair}_2.fq.gz \
        -O ../filtered_reads/${SAMPLE}/Filtered_${ReadPair}_2.fq.gz \
        --trim_front1 10 \
        --trim_front2 10 \
        --detect_adapter_for_pe \
        --dedup
        
        #Rename / move fastp reports to fastp_log directory
        mv fastp.html ../fastp_logs/${ReadPair}.html
        mv fastp.json ../fastp_logs/${ReadPair}.json

        #Send fastp reports to s3 bucket for storage
        ~/aws-cli/bin/aws s3 cp ../fastp_logs/${ReadPair}.html s3://norfolkhybrids/fastp_reports/
        ~/aws-cli/bin/aws s3 cp ../fastp_logs/${ReadPair}.json s3://norfolkhybrids/fastp_reports/

        #Send filtered reads to s3 bucket
        ~/aws-cli/bin/aws s3 cp ../filtered_reads/${SAMPLE}/Filtered_${ReadPair}_1.fq.gz s3://norfolkhybrids/filtered_reads/${SAMPLE}/
        ~/aws-cli/bin/aws s3 cp ../filtered_reads/${SAMPLE}/Filtered_${ReadPair}_2.fq.gz s3://norfolkhybrids/filtered_reads/${SAMPLE}/

    done

    #Move out of temp directory and delete it
    cd ../
    rm -r temp

    #Add entry to log file
    echo "read filtering of sample" $SAMPLE "complete!" >> filtered_reads/read_filtering.log

done
