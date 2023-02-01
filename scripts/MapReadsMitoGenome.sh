mkdir mapping/Zlat_mitogenome
#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   ADD
# Author:   Ash Sendell-Price
# Date:     18/01/2023
# System:   vettel (warwick)
# Summary:  Maps filtered novogene reads to Zosterops lateralis mitogenome assembly, merges bams,
#           calculates sequencing depth at each position.
#################################################################################################

#Load conda environment
mamba activate NorfolkIslandHybridisation

#Create directories to store sample bam files
mkdir mapping
mkdir mapping/Zlat_mitogenome
mkdir mitogenome_fastas

#Index assembly using bwa-mem2
bwa-mem2 index assemblies/NC_029146.1_ZLat_mitogenome.fasta.gz

#For each of the novogene sequenced samples in resources/novogene_samples.txt
#do the following:
for SAMPLE in $(cat resources/novogene_samples.txt)
do
    #Make a temporary directory and move into it
    mkdir mapping_mito_temp
    cd mapping_mito_temp

    #For each read pair do the following ...
    for ReadPair in $(ls ../novogene_filtered_reads/${SAMPLE}/Filtered_${SAMPLE}_*_1.fq.gz | rev | cut -d "_" -f2- | rev)
    do
        BASE_NAME=$(basename $ReadPair)
        bwa-mem2 mem \
        ../assemblies/NC_029146.1_ZLat_mitogenome.fasta.gz \
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
    samtools sort ${SAMPLE}.bam > ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam

    #Index bam file
    samtools index ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam

    #Calculate depth at each site (incl. missing)
    samtools depth -aa ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam > ../mapping/Zlat_mitogenome/${SAMPLE}.depth

    #Output consensus mitogenome sequence for sample
    ../bin/angsd/angsd -i ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam \
    -minQ 20 -minMapQ 20 -uniqueOnly -setMinDepth 10 \
    -doCounts 1 -dumpCounts 1 -doFasta 4 \
    -out mitogenome_fastas/${SAMPLE}

    #Transfer full bam and index to s3 bucket
    ~/aws-cli/bin/aws s3 cp ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam s3://norfolkhybrids/mapping/Zlat_mitogenome/
    ~/aws-cli/bin/aws s3 cp ../mapping/Zlat_mitogenome/${SAMPLE}.sorted.bam.bai s3://norfolkhybrids/mapping/Zlat_mitogenome/
    ~/aws-cli/bin/aws s3 cp ../mapping/Zlat_mitogenome/${SAMPLE}.depth s3://norfolkhybrids/mapping/Zlat_mitogenome/

    #Tidy up
    cd ../
    rm -r mapping_mito_temp

done
