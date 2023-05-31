#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   NGS_Admix.sh
# Author:   Ash Sendell-Price
# Date:     02/05/2023
# System:   vettel (warwick)
# Summary:  ADD
#################################################################################################

# Get current working directory "NorfolkIslandHybridization"
PROJECT_DIRECTORY=$(pwd)

# Make directory for analysis and move into it
mkdir ${PROJECT_DIRECTORY}/admixture
cd ${PROJECT_DIRECTORY}/admixture

# Remove outgroup samples from beagle file
zcat ../genotype_likelihoods_ANGSD/NorfolkIsland_NewZealand_plus_Zvirens_minMAF0.05.beagle.gz \
| cut -f 1-105 | bgzip > NorfolkIsland_NewZealand_Zlat_Zten_only.beagle.gz

#Run NGSadmix testing K range 2 to 5, performing 5 replicates per K
for K in $(seq 2 2)
do
    for RUN in $(seq 1 5)
    do
        ${PROJECT_DIRECTORY}/bin/NGSadmix -likes NorfolkIsland_NewZealand_Zlat_Zten_only.beagle.gz \
        -K $K -minMaf 0.05 -o NorfolkIsland_NewZealand_Zlat_Zten_only_K${K}_run${RUN}
    done
done

