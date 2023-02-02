library(tidyverse)
library(adegenet) 
library(ade4) 

gl = load("./gl.RData") 

pop.list <- c("SUD", "SUD", "SUD", "SUD", "SUD", 
              "VLA", "VLA", "VLA", "VLA", "VLA", 
              "TON", "TON", "TON", "TON", "TON", 
              "SUD","SUD","SUD","SUD", 
              "VLA",
              "TON", 
              "CHA", "CHA", "CHA", "CHA", "CHA", 
              "KES", "KES", "KES", 
              "CKM", "CKM", "CKM", "CKM", "CKM", 
              "KES", 
              "CKM","CKM","CKM")
pop(gl) <- pop.list  

pca1 <- glPca(gl, 
              parallel = TRUE) 

save(pca1, file = "./pca1.RData")