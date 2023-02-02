library(vcfR)                                                                    
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
