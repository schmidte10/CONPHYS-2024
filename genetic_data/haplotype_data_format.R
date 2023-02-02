#--- lirbary ---# 
library(vcfR) 
library(geneHapR)
#--- script ---# 
snps <- read.vcfR("./genetic_vcf/SCH11429_finalgatk.vcf") 

haplotype_data <-  geneHapR::vcf2hap(snps, 
                                     hapPrefix = "H", 
                                     filter_Chr = FALSE, 
                                     filter_POS = FALSE, 
                                     hetero_remove = TRUE, 
                                     na_drop =TRUE) 


save(haplotype_data, file = "./haplotype_data.RData")
save(haplotype_data, file = "./haplotype_data.txt")