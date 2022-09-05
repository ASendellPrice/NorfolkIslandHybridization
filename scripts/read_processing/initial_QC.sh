#!/bin/bash

#####################################################################
# Parameters
#####################################################################

#Set path to projectc directory
PROJECT=/data/Users/Ash/NorfolkIslandHybridization/

#Set path to file specifying read pair info
READ_PAIR_INFO=${PROJECT}/scripts/read_processing/read_info.txt


#####################################################################
#####################################################################


#If directory "mapping" does not exist then create it
if [[ ! -d ${PROJECT}/data/raw_reads/read_qc ]]
then
    mkdir ${PROJECT}/data/raw_reads/read_qc
fi

#Move into mapping directory
cd ${PROJECT}/data/raw_reads/read_qc


#For each line in READ_PAIR_INFO do the following
for READ_PAIR in $(cat $READ_PAIR_INFO)
do
    SAMPLE_ID=$(echo $READ_PAIR | cut -f 1)
done



fastqc
