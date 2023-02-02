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

module load R
singularity run $SING/R-4.1.2.sif R

Rscript -e "library(vcfR)"                                                                    
library(stringr)
library(readr) 
library(dartR)

snps <- read.vcfR("./genetic_vcf/SCH11429_finalgatk.vcf")                        
snps_gl <- vcfR2genlight(snps) 

pops <- snps_gl@ind.names %>% str_remove(pattern="_.*") 
inds <- snps_gl@ind.names %>% str_remove(pattern = ".*_") %>% parse_number() 

snps_gl@pop <- as.factor(pops) 
snps_gl@ind.names <- as.chacaracter(inds) 

# re-arrange the individuals 
snps_gl_reorder <- snps_gl[order(as.numeric(snps_gl@ind.names))]
snps_gl_reorder@ind.names

gl <- snps_gl_reorder
#We do a compliance check so the metadata slot is filled
gl <- gl.compliance.check(gl)
save(gl, file = "./gl.RData")

