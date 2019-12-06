setwd("~/Desktop/BranchyData/data/")
library(tidyverse)
options(warn = -1)

rm(list=ls())
#################the QTL effect for each site 
## we need to run the scan1coef for each trait at each site, need to reformat the data file
###since running single site, we need to build the yaml file for single site

#pheno for each CO2 level
phenos = read.csv('Phenotypes_Cleaned.csv')

pheno_single = vector("list",length = 0)
z=0
for (tr in unique(phenos$Treatment)){
  z=z+1
  tmp = phenos %>% filter(Treatment==tr) %>% dplyr::select(-Treatment)
  colnames(tmp)= c('id', paste0(colnames(phenos)[-c(1:2)], "_",tr))
  pheno_single[[z]] = tmp
}

pheno = pheno_single %>% bind_cols() %>% dplyr::select(-matches("id[0-9]"))

##pheno for differences between high/low CO2 and the control
pheno = pheno %>% mutate(AmountC_high = AmountC_elevated - AmountC_ambient, AmountC_low = AmountC_subamb - AmountC_ambient,
                              AmountN_high = AmountN_elevated - AmountN_ambient, AmountN_low = AmountN_subamb - AmountN_ambient,
                              CN_high = CN_elevated - CN_ambient, CN_low = CN_subamb - CN_ambient,
                              d13C_high = d13C_elevated - d13C_ambient, d13C_low = d13C_subamb - d13C_ambient,
                              d15N_high = d15N_elevated - d15N_ambient, d15N_low = d15N_subamb - d15N_ambient)

## pheno for averages across CO2 level
pheno_avg = phenos %>% dplyr::select(-Treatment) %>% group_by(genotype) %>% 
  summarise_all(mean, na.rm=TRUE) %>% dplyr::rename(id = genotype)

pheno_all = pheno %>% left_join(pheno_avg) 

write.csv(pheno_all,"RIL_Single_phe.csv",row.names = F,quote = F)

#geno 
n= length(unique(phenos$Treatment))
genos = read.csv('RIL_Multiple_gen.csv')
geno = genos[1:(nrow(genos)/n),]
geno = geno %>% mutate(id = pheno$id)
write.csv(geno,"RIL_Single_gen.csv",row.names = F,quote = F)

#map 
file.copy("RIL_Multiple_gmap.csv","RIL_Single_gmap.csv")

#yaml (control) file 
cat('# Data for 165 RIL_Single in Single Site:\n',file = 'RIL_Single.yaml')
cat('crosstype: riself\n',file = 'RIL_Single.yaml',append = T)
cat('geno: RIL_Single_gen.csv\n',file = 'RIL_Single.yaml',append = T)
cat('pheno: RIL_Single_phe.csv\n',file = 'RIL_Single.yaml',append = T)
cat('gmap: RIL_Single_gmap.csv\n',file = 'RIL_Single.yaml',append = T)
cat('alleles:\n- A\n- B\n',file = 'RIL_Single.yaml',append = T)
cat('genotypes:\n',
    'A: 1\n',
    'B: 2\n',file = 'RIL_Single.yaml',append = T)
cat('na.strings:\n',
    "- '-'\n",
    '- NA\n',file = 'RIL_Single.yaml',append = T)






