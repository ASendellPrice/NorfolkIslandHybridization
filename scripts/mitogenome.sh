#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   ADD
# Author:   Ash Sendell-Price
# Date:     18/01/2023
# System:   vettel (warwick)
# Summary:  ADD
#################################################################################################

#Load conda environment
mamba activate NorfolkIslandHybridisation

#Create directories to store sample bam files
mkdir mitogenomes

#For each of the novogene sequenced samples in resources/novogene_samples.txt
#do the following:
for SAMPLE in $(cat resources/novogene_samples.txt)
do
    #Generate mitochonrial fasta sequence using angsd - output to temp file
    bin/angsd/angsd -i mapping/Zlat_2_Tgutt_ZW/${SAMPLE}.sorted.bam \
    -minQ 20 -minMapQ 20 -uniqueOnly -setMinDepth 5 \
    -doCounts 1 -dumpCounts 1 -doFasta 4 \
    -r PseudoCM016613.2_Taeniopygia_guttata_isolate_Black17_mitochondrion \
    -out mitogenomes/${SAMPLE}
done
