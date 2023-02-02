#!/bin/bash
#PBS -j oe
#PBS -N APOLYHPC
#PBS -l select=1:ncpus=32:mem=100gb
#PBS -l walltime=18:00:00
#PBS -m ae
#PBS -M elliott.schmidt@my.jcu.edu.au

echo "------------------------------------------------------"
echo "PBS: Submitted to $PBS_QUEUE@$PBS_O_HOST"
echo "PBS: Working directory is $PBS_O_WORKDIR"
echo "PBS: Job identifier is $PBS_JOBID"
echo "PBS: Job name is $PBS_JOBNAME$no"
echo "------------------------------------------------------"

module load anaconda3 
source $CONDA_PROF/conda.sh 
conda activate pear-0.9.6 

pear -f ./sequencing_data/CSUD006_S1_L004_R1_001.fastq.gz -r ./sequencing_data/CSUD006_S1_L004_R2_001 -o ./sequencing_data/merge/CSUD_006 
pear -f ./sequencing_data/CSUD008_S1_L004_R1_001.fastq.gz -r ./sequencing_data/CSUD008_S1_L004_R2_001 -o ./sequencing_data/merge/CSUD_008 