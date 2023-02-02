
library(vcfR)                                                                    # load 'vcfR' package
library(stringr)
library(readr) 
library(dartR)

snps <- read.vcfR("./genetic_vcf/SCH11429_finalgatk.vcf")                        # read vcf file 
snps_gl <- vcfR2genlight(snps) 

pops <- snps_gl@ind.names %>% str_remove(pattern="_.*") 
inds <- snps_gl@ind.names %>% str_remove(pattern = ".*_") %>% parse_number() 

snps_gl@pops <- as.factor(pops) 
snps_gl@ind.names <- as.character(inds) 

# re-arrange the individuals 
snps_gl_reorder <- snps_gl[order(as.numeric(snps_gl@ind.names))] 
snps_gl_reorder@ind.names 
gl <- snps_gl_reorder 

gl <- gl.compliance.check(gl)

save(gl, file = "./gl.RData")