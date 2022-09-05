#!/bin/bash

#####################################################################
# Parameters
#####################################################################

#Set path to projectc directory
PROJECT=/data/Users/Ash/NorfolkIslandHybridization/

#Set path to file specifying read pair info
READ_PAIR_INFO=${PROJECT}/


#####################################################################
#####################################################################


#If directory "mapping" does not exist then create it
if [[ ! -d ${PROJECT}/read_qc ]]
then
    mkdir ${PROJECT}/read_qc
fi

#Move into mapping directory
cd ${PROJECT}/read_qc




fastqc
