#!/bin/bash
#PBS -j oe
#PBS -N APOLYHPC
#PBS -l select=1:ncpus=32:mem=16gb
#PBS -l walltime=18:00:00
#PBS -m ae
#PBS -M elliott.schmidt@jcu.edu.au

cd $PBS_O_WORKDIR

echo "------------------------------------------------------"
echo "PBS: Submitted to $PBS_QUEUE@$PBS_O_HOST"
echo "PBS: Working directory is $PBS_O_WORKDIR"
echo "PBS: Job identifier is $PBS_JOBID"
echo "PBS: Job name is $PBS_JOBNAME$no"
echo "------------------------------------------------------"

module load singularity
singularity run $SING/R-4.1.2.sif R --vanilla
R -f~/genetic_vcf/local_adapt.R

