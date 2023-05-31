#!/bin/bash

#################################################################################################
# Project:  Norfolk Island hybridization study
# Script:   Calc_Dstats_Dstuite_GLs.sh
# Author:   Ash Sendell-Price
# Date:     11/04/2023
# System:   vettel (warwick)
# Summary:  ADD
#################################################################################################

# Get current working directory "NorfolkIslandHybridization"
PROJECT_DIRECTORY=$(pwd)

# Make directory for analysis and move into it
mkdir ${PROJECT_DIRECTORY}/Dsuite_introgression
cd ${PROJECT_DIRECTORY}/Dsuite_introgression

#Set population assignment file.
#This file has two columns: sampleID and population
POP_ASSIGNMENT=SETS_Zvirens.txt

#Set trios file
TRIOS=test_trios.txt

#Run dtrios to calcualte chromosome-wide D statistics
${PROJECT_DIRECTORY}/bin/Dsuite/Build/Dsuite Dtrios \
${PROJECT_DIRECTORY}/genotype_likelihoods_ANGSD/NorfolkIsland_NewZealand_plus_Zvirens_minMAF0.05_SampleIDsUpdated.vcf.gz \
$POP_ASSIGNMENT --use-genotype-probabilities

#Run dinvestigate to calculate windowed d-statistics across the chromosome
${PROJECT_DIRECTORY}/bin/Dsuite/Build/Dsuite Dinvestigate \
${PROJECT_DIRECTORY}/genotype_likelihoods_ANGSD/NorfolkIsland_NewZealand_plus_Zvirens_minMAF0.05_SampleIDsUpdated.vcf.gz \
$POP_ASSIGNMENT $TRIOS --use-genotype-probabilities --window=500,100 --run-name OutgroupZvirens

