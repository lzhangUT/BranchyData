setwd("~/Desktop/BranchyData/0-RawData/")
library(tidyverse)
library(corrplot)
options(warn = -1)

rm(list = ls())
###read in the biomass data and get the average for each environment each ril

Bio_raw = read.csv('Brachy RIL x CO2 biomass 27 April 2016.csv')

Bio = Bio_raw %>% filter(!notes =="corn") %>% select(Trt, RIL, matches("mass|ht")) %>% 
  rowwise() %>% mutate(biomass = sum(leaf_mass, stem_mass, seed_mass)) %>% filter(!RIL=="blank")

Bio_Avg = Bio %>% group_by(RIL, Trt)%>% summarize_all(list(~mean(., na.rm = T))) %>% 
  ungroup() %>% as.data.frame() %>% mutate_if(is.numeric, round, 2) %>% mutate_at(vars(RIL), as.character) %>%
  mutate(RIL= ifelse(nchar(RIL)==3, paste0("RIL",RIL), RIL)) %>%
  mutate(RIL= ifelse(nchar(RIL)==2, paste0("RIL0",RIL), RIL)) %>%
  mutate(RIL= ifelse(nchar(RIL)==1, paste0("RIL00",RIL), RIL)) %>%
  dplyr::rename(Treatment = Trt) %>% mutate(RIL=str_replace(RIL,"33333","bd3-1")) %>% 
  mutate(RIL=str_replace(RIL,"212121","bd21"))

##read in the raw CN data
CN_raw = readxl::read_excel('Juenger BRPRIL_AUS_2015 0317.xls', skip = 6)
SampleID = readxl::read_excel('BRPRIL_2015_Isotope_Sample_Pooling_Processing_Order 2016-11-07.xlsx')

SampleID = SampleID %>% dplyr::select(Tissue, `RIL Line`, Treatment,`Plate&Postion`) %>% 
  dplyr::rename(RIL = `RIL Line`) %>% filter(Tissue =="leaves") %>% 
  separate(Treatment, c("Treatment", 'other'))

CN = CN_raw %>% filter(str_detect(`Sample ID`,"leaves")) %>% 
  dplyr::select("Sample ID", "d13C","d15N","C Amount (ug)","N Amount (ug)")%>%
  separate("Sample ID",c("Plate","Position","Tissue","RIL","other")) %>% 
  unite("Plate&Postion",c("Plate","Position"), sep = '-') %>%
  mutate(RIL, RIL=ifelse(is.na(other), RIL,paste(RIL,other, sep = "-"))) %>% 
  dplyr::select(-other) %>% mutate(CN = `C Amount (ug)`/`N Amount (ug)`) %>%
  dplyr::rename(AmountC =`C Amount (ug)`, AmountN = `N Amount (ug)` )

CN = CN %>% left_join(SampleID)%>% filter(complete.cases(.)) %>%
  mutate(RIL, RIL= ifelse(nchar(RIL)==3, paste0("RIL",RIL), RIL)) %>%
  mutate(RIL, RIL= ifelse(nchar(RIL)==2, paste0("RIL0",RIL), RIL)) %>%
  mutate(RIL, RIL= ifelse(nchar(RIL)==1, paste0("RIL00",RIL), RIL)) %>% 
  dplyr::select(-other, -`Plate&Postion`, -Tissue)

###merge CN with BIO data
phenotypes_raw = Bio_Avg %>% left_join(CN) %>% dplyr::rename(genotype = RIL) 
phenotypes_log = phenotypes_raw %>% dplyr::select(genotype, Treatment, matches("mass"))%>% 
  mutate_at(vars(matches("mass")), log) %>% rename(leaf_mass_log = leaf_mass, stem_mass_log = stem_mass,
                                                   seed_mass_log = seed_mass, biomass_log = biomass)

phenotypes_raw_log = phenotypes_raw %>% left_join(phenotypes_log) %>%
  select(genotype, Treatment,matches("mass"),matches("ht"), d13C:CN)

##############################################################
#####plot histogram and Deal with outlier for phenotype file##
##############################################################
phenos = phenotypes_raw_log %>% gather(Trait, Value, -genotype, -Treatment)
phenos$Treatment = factor(phenos$Treatment, levels = c('elevated','ambient','subamb'))

pdf('../3-Figures/0-Histogram with Raw Data.pdf', onefile = T, width = 10)

for (tr in unique(phenos$Treatment)){
  par(mfrow=c(4,4))
  for (phe in unique(phenos$Trait)){
    tmp = phenos %>% filter(Trait==phe & Treatment==tr)
    hist(tmp$Value, main = tr, xlab = phe)
  }
}

###plot in ggplot2
ggplot(phenos, aes(Value)) + geom_histogram(bins = 30)+facet_grid(Treatment~Trait, scales = "free")+ theme_bw()
dev.off()

###Outliers need to be removed for each trait
pdf('../3-Figures/0-Histogram with Outlier Removed.pdf', onefile = T, width = 10)

phenos_clean = NULL
for (tr in unique(phenos$Treatment)){
  par(mfrow=c(4,4))
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


phenos_clean2 = phenos_clean %>% spread(Trait,Value) %>% arrange(Treatment)%>% 
  select(genotype, Treatment, matches("mass"),matches("ht"),  matches("Amount"),d13C, d15N, CN) %>%
  mutate(genotype=str_replace(genotype,"bd3-1","Bd_3_1")) %>% 
  mutate(genotype=str_replace(genotype,"bd21","Bd_21"))
  
write.csv(phenos_clean2, "../1-ProcessedData/Phenotypes_Cleaned.csv", row.names = F)

######phenotypic correlation between traits
pdf("../3-Figures/0-PhenotypicCorrelation.pdf", width = 10, height = 8)
for (tr in unique(phenos_clean2$Treatment)){
  tmp = phenos_clean2 %>% select(Treatment, leaf_mass, stem_mass, seed_mass, biomass,ht_flag,ht_flower, AmountC:CN)%>% 
    filter(Treatment==tr)%>% select(-Treatment) %>% filter(complete.cases(.))
  M = cor(tmp)
  corrplot.mixed(M, title=tr,mar = c(1,1,2.1,1))
  
}
phenos_cor = phenos_clean2 %>% select(leaf_mass, stem_mass, seed_mass, biomass, AmountC:CN) %>% 
  filter(complete.cases(.))
M = cor(phenos_cor)
corrplot.mixed(M, title="All Treatments",mar = c(1,1,2.1,1))

dev.off()



