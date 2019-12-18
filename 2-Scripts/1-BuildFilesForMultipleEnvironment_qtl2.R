setwd("~/Desktop/BranchyData/1-ProcessedData/")
library(tidyverse)
options(warn = -1)

rm(list = ls())

#To run in rqtl2 with covariates (in this case, different CO2 levels), the following files are needed:
#   geno: RIL_Multiple_gen.csv
#   pheno: RIL_Multiple_phe.csv
#   covar: RIL_Multiple_cov.csv 
#   gmap: RIL_Multiple_gmap.csv
#   ymal: RIL_Multiple.yaml (# the control file)


###. working on the genotype and map file 
cross = read.csv("../0-RawData//Brachypodium_RQTL data.csv", na.strings = c("","NA"))
##remove the original phenotype data there
cross2 = cross %>% dplyr::select(-matches("amb|elev|sub") ) %>% filter_all(any_vars(!is.na(.)))

id = cross2 %>% dplyr::select(genotype) %>% filter(!is.na(genotype)) 

###line up the phentoype file with the genotype files
phenos = read.csv("Phenotypes_Cleaned.csv")
n = length(unique(phenos$Treatment))

phenos_cov = NULL
for (i in unique(phenos$Treatment)){
  tmp = id %>% left_join(phenos %>% filter(Treatment==i))%>% mutate(Treatment=i)
  phenos_cov = rbind(tmp,phenos_cov)
}
write.csv(phenos_cov, "Phenotypes_LinedUpWithGenotypes.csv", row.names = F)

##pheno file
pheno = phenos_cov %>% mutate(id=1:n()) %>% dplyr:: select(id,biomass:CN)
write.csv(pheno,"RIL_Multiple_phe.csv",row.names = F,quote = F)

###covarite file
covar = phenos_cov %>% mutate(id=1:n()) %>% dplyr:: select(id, Treatment)%>% 
  mutate_at(vars(Treatment), as.factor) %>% mutate_at(vars(Treatment), as.numeric)
write.csv(covar,"RIL_Multiple_covar.csv",row.names = F,quote = F)

#geno, repeat the geno file n times, n is the number of sites
tmp_geno = cross2[-c(1:2),-1]
geno = bind_rows(replicate(n,tmp_geno, simplify = FALSE)) %>% mutate(id= 1:n()) %>% dplyr::select(id, colnames(tmp_geno))
write.csv(geno,"RIL_Multiple_gen.csv",row.names = F,quote = F)

##map
map_raw = cross2 %>% select(-genotype) %>% slice(1:2)
map = as.data.frame(t(map_raw))
map = map %>% rownames_to_column('marker') %>% rename(chr = "V1", pos="V2")
write.csv(map,"RIL_Multiple_gmap.csv",row.names = F,quote = F)

#yaml (control) file
cat('# Isotope and Biomass data for 165 RIL in Multiple Treatment:\n',file = 'RIL_Multiple.yaml')
cat('crosstype: riself\n',file = 'RIL_Multiple.yaml',append = T)
cat('geno: RIL_Multiple_gen.csv\n',file = 'RIL_Multiple.yaml',append = T)
cat('pheno: RIL_Multiple_phe.csv\n',file = 'RIL_Multiple.yaml',append = T)
cat('covar: RIL_Multiple_covar.csv\n',file = 'RIL_Multiple.yaml',append = T)
cat('gmap: RIL_Multiple_gmap.csv\n',file = 'RIL_Multiple.yaml',append = T)
cat('alleles:\n- A\n- B\n',file = 'RIL_Multiple.yaml',append = T)
cat('genotypes:\n',
    'A: 1\n',
    'B: 2\n',file = 'RIL_Multiple.yaml',append = T)
cat('Treatment:\n',
    'covar: Treatment\n',
    'ambient: 1\n',
    'elevated: 2\n',
    'subamb: 3\n', file = 'RIL_Multiple.yaml',append = T)
cat('na.strings:\n',
    "- '-'\n",
    '- NA\n',file = 'RIL_Multiple.yaml',append = T)


#############END


