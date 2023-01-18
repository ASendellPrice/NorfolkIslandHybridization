#!/bin/bash

#Install following conda packages
conda install -c bioconda samtools

#Sync assemblies directory
aws s3 sync s3://norfolkhybrids/assemblies ./assemblies

#Copy lists of samples names and chrom names from s3 bucket
aws s3 cp s3://norfolkhybrids/sample_ids.txt ./
aws s3 cp s3://norfolkhybrids/chrom.names.txt ./

#Download bwa-mem2 executable
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
| tar jxf -

#Index psuedoassembly using bwa-mem2
REF=/home/ec2-user/assemblies/Zlat_2_Tgutt_ZW.fasta.gz
BWA=/home/ec2-user/bwa-mem2-2.0pre2_x64-linux/bwa-mem2
$BWA index $REF
aws s3 sync ./assemblies s3://norfolkhybrids/assemblies

#For each sample in list do the following ... 
for SAMPLE in $(cat sample_ids.txt | tail -n +2)
do
    #Make a directory for sample and move into it
    mkdir $SAMPLE
    cd $SAMPLE

    #Copy over samples' filtered reads from s3 bucket and save them in a directory called "filtered_reads"
    aws s3 sync s3://norfolkhybrids/filtered_reads/${SAMPLE} ./filtered_reads

    #Make directory for bam files
    mkdir BAMs




##### USE BWA MEM2 on vettel bwa-mem2


    #For each read pair do the following ...
    for ReadPair in `ls filtered_reads/Filtered_${SAMPLE}_*_1.fq.gz | cut -f1,2,3,4,5,6 -d'_'`
    do
        BASE_NAME=$(basename $ReadPair)
        $BWA mem -t 30 $REF \
        ${ReadPair}_1.fq.gz \
        ${ReadPair}_2.fq.gz \
        -R "@RG\tID:${BASE_NAME}\tSM:${SAMPLE}" \
        | samtools view -bS - \
        | samtools sort - > BAMs/${BASE_NAME}.bam
    done

    #If number of BAMs is greater than 1 then merge bams into a single file else
    #just rename the bam file
    if [ $(ls BAMs/Filtered_${SAMPLE}*.bam | wc -l) -gt 1 ]
    then
        samtools merge BAMs/${SAMPLE}.bam BAMs/Filtered_${SAMPLE}*.bam -f
        rm BAMs/Filtered_${SAMPLE}*.bam
    else
        mv BAMs/Filtered_${SAMPLE}*.bam BAMs/${SAMPLE}.bam
    fi

    #Sort bam file and remove unsorted version
    samtools sort BAMs/${SAMPLE}.bam > BAMs/${SAMPLE}.sorted.bam
    rm BAMs/${SAMPLE}.bam

    #Index bam file
    samtools index BAMs/${SAMPLE}.sorted.bam

    #For each psuedochromosome do the following ...
    for CHROM in $(cat /home/ec2-user/chrom.names.txt)
    do
        #Extract read within focal psuedochromosome
        samtools view -b BAMs/${SAMPLE}.sorted.bam ${CHROM} > BAMs/${CHROM}.${SAMPLE}.sorted.bam
        #Index the file
        samtools index BAMs/${CHROM}.${SAMPLE}.sorted.bam
        #Calculate depth
        samtools depth BAMs/${CHROM}.${SAMPLE}.sorted.bam > BAMs/${CHROM}.${SAMPLE}.depth.txt
        #Transfer files to s3 bucket
        aws s3 cp BAMs/${CHROM}.${SAMPLE}.sorted.bam s3://norfolkhybrids/chrom_bams/${CHROM}/
        aws s3 cp BAMs/${CHROM}.${SAMPLE}.sorted.bam.bai s3://norfolkhybrids/chrom_bams/${CHROM}/
        aws s3 cp BAMs/${CHROM}.${SAMPLE}.depth.txt s3://norfolkhybrids/chrom_bams/${CHROM}/
    done

    #Transfer full bam and index to s3 bucket
    aws s3 cp BAMs/${SAMPLE}.sorted.bam s3://norfolkhybrids/sample_bams/
    aws s3 cp BAMs/${SAMPLE}.sorted.bam.bai s3://norfolkhybrids/sample_bams/

    #Tidy up
    cd ../
    rm -r $SAMPLE
done

