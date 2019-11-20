setwd("~/Desktop/BranchyData/data/")
library(tidyverse)
options(warn = -1)

rm(list = ls())
###1. working on formating phenotype file
CN_raw = readxl::read_excel('Juenger BRPRIL_AUS_2015 0317.xls', skip = 6)

CN = CN_raw %>% filter(str_detect(`Sample ID`,"leaves")) %>% 
  dplyr::select("Sample ID", "d13C","d15N","C Amount (ug)","N Amount (ug)")%>%
  separate("Sample ID",c("Plate","Position","Tissue","RIL","other")) %>% 
  unite("Plate&Postion",c("Plate","Position"), sep = '-') %>%
  mutate(RIL, RIL=ifelse(is.na(other), RIL,paste(RIL,other, sep = "-"))) %>% 
  dplyr::select(-other) %>% mutate(CN = `C Amount (ug)`/`N Amount (ug)`) %>%
  dplyr::rename(AmountC =`C Amount (ug)`, AmountN = `N Amount (ug)` )

##combining with the environmental data 
SampleID = readxl::read_excel('BRPRIL_2015_Isotope_Sample_Pooling_Processing_Order 2016-11-07.xlsx')

tmp_pheno = CN %>% left_join( SampleID %>% dplyr::select(Tissue, `RIL Line`, Treatment,`Plate&Postion`) %>% 
                            dplyr::rename(RIL = `RIL Line`) %>% filter(Tissue =="leaves") %>% 
            separate(Treatment, c("Treatment", 'other')) ) %>% filter(complete.cases(.)) %>%
            mutate(RIL, RIL= ifelse(nchar(RIL)==3, paste0("RIL",RIL), RIL)) %>%
            mutate(RIL, RIL= ifelse(nchar(RIL)==2, paste0("RIL0",RIL), RIL)) %>%
            mutate(RIL, RIL= ifelse(nchar(RIL)==1, paste0("RIL00",RIL), RIL)) %>% 
           dplyr::select(-other, -`Plate&Postion`, -Tissue) ###this phenotype file needs to be lined up with genotype file

###2. working on the genotype and map file 
cross = read.csv("Brachypodium_RQTL data.csv", na.strings = c("","NA"))
##remove the original phenotype data there
cross2 = cross %>% dplyr::select(-matches("amb|elev|sub") ) %>% filter_all(any_vars(!is.na(.))) %>% 
  filter(genotype != "Bd_3_1") %>% filter(genotype!= "Bd_21")

id = cross2 %>% dplyr::select(genotype) %>% filter(!is.na(genotype)) 

###line up the phentoype file with the genotype files
n = length(unique(tmp_pheno$Treatment))

phenos_cov = NULL
for (i in unique(tmp_pheno$Treatment)){
  tmp = id %>% left_join(tmp_pheno %>% filter(Treatment==i), by=c("genotype"="RIL"))
  phenos_cov = rbind(tmp,phenos_cov)
}
phenos_cov = phenos_cov %>% fill(Treatment)

##############################################################
#####plot histogram and Deal with outlier for phenotype file##
##############################################################
phenos = phenos_cov  %>% gather(Trait, Value, -genotype, -Treatment)
phenos$Treatment = factor(phenos$Treatment, levels = c('elevated','ambient','subamb'))

pdf('../figures/Histogram with Raw Data.pdf', onefile = T)
par(mfrow=c(3,5))
for (tr in unique(phenos$Treatment)){
  for (phe in unique(phenos$Trait)){
    tmp = phenos %>% filter(Trait==phe & Treatment==tr)
    hist(tmp$Value, main = tr, xlab = phe)
  }
}

###plot in ggplot2
ggplot(phenos, aes(Value)) + geom_histogram(bins = 30)+facet_grid(Treatment~Trait, scales = "free")+ theme_bw()
dev.off()

###Outliers need to be removed for each trait
pdf('../figures/Histogram with Outlier Removed.pdf', onefile = T)
par(mfrow=c(3,5))
phenos_clean = NULL
for (tr in unique(phenos$Treatment)){
  for (phe in unique(phenos$Trait)){
  tmp = phenos %>% filter(Trait==phe & Treatment==tr)
  outliers = boxplot.stats(tmp$Value,coef = 1.5)$out
  tmp$Value[tmp$Value %in% outliers] = NA
  hist(tmp$Value, main = tr, xlab = phe)
  phenos_clean = rbind(phenos_clean,tmp)
  }
}

ggplot(phenos_clean, aes(Value)) + geom_histogram(bins = 30)+facet_grid(Treatment~Trait, scales = "free")+ theme_bw()
dev.off()

phenos_clean2 = phenos_clean %>% spread(Trait,Value) %>% arrange(Treatment)
write.csv(phenos_clean2, "Phenotypes_Cleaned.csv", row.names = F)
###outlier removed, now working to build the files needed for QTL mapping

#To run in rqtl2 with covariates (in this case, different CO2 levels), the following files are needed:
#   geno: RIL_Multiple_gen.csv
#   pheno: RIL_Multiple_phe.csv
#   covar: RIL_Multiple_cov.csv 
#   gmap: RIL_Multiple_gmap.csv
#   ymal: RIL_Multiple.yaml (# the control file)


##pheno file
pheno = phenos_clean2 %>% mutate(id=1:n()) %>% dplyr:: select(id, "d13C","d15N","AmountC","AmountN","CN")
write.csv(pheno,"RIL_Multiple_phe.csv",row.names = F,quote = F)

###covarite file
covar = phenos_cov %>% mutate(id=1:n()) %>% dplyr:: select(id, Treatment)%>% 
  mutate_at(vars(Treatment), as.factor) %>% mutate_at(vars(Treatment), as.numeric)
write.csv(covar,"RIL_Multiple_covar.csv",row.names = F,quote = F)

#geno, repeat the geno file n times, n is the number of sites
tmp_geno = cross2[,-1]
geno = bind_rows(replicate(n,tmp_geno, simplify = FALSE)) %>% mutate(id= 1:n()) %>% dplyr::select(id, colnames(tmp_geno))
write.csv(geno,"RIL_Multiple_gen.csv",row.names = F,quote = F)

##map
map_raw = cross %>% dplyr::select(-matches("amb|elev|sub") ) %>% filter_all(any_vars(!is.na(.))) %>% 
  select(-genotype) %>% slice(1:2)
map = as.data.frame(t(map_raw))
map = map %>% rownames_to_column('marker') %>% rename(chr = "V1", pos="V2")
write.csv(map,"RIL_Multiple_gmap.csv",row.names = F,quote = F)

#yaml (control) file
cat('# Isotope Data for about 165 RIL in Multiple Treatment:\n',file = 'RIL_Multiple.yaml')
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
