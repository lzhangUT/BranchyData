setwd("~/Desktop/BranchyData/data/")
library(qtl)
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(tidyverse)
options(warn = -1)

rm(list=ls())

###using Rqtl2 for RIL population
cross2 = read_cross2("RIL_Single.yaml")
map = insert_pseudomarkers(cross2$gmap, step=1)
pr_single = calc_genoprob(cross2, map, error_prob=0.002,cores = 2)
kinship_loco_single = calc_kinship(pr_single, "loco", cores=2)
####
pheno = cross2$pheno
pheid = colnames(pheno)[-1]
EffectMatrix = vector('list', length=0)
t = 0
#######for each trait run all the chromosome for the effect
for (phe in pheid){
  t=t+1
  z=0
  eff_out = vector('list',length = 0)
  se_out = vector('list',length = 0)
  for (chr in chrnames(cross2)){
    z = z+1
    eff = scan1coef(pr_single[,chr],pheno[,phe],kinship = kinship_loco_single[chr], se=T)
    se = attr(eff,'SE')
    eff_out[[z]] =  eff %>% as.data.frame() %>% rownames_to_column('MARKER') 
    se_out[[z]] =  se %>% as.data.frame() %>% rownames_to_column('MARKER')
  }
  eff_out = bind_rows(eff_out) 
  se_out = bind_rows(se_out)
  
  se_out2 = se_out %>% dplyr::select(-intercept) %>% gather('CROSS','SE',-MARKER) 
  eff_out2 = eff_out %>% dplyr::select(-intercept) %>% gather('CROSS','EFF',-MARKER) %>% mutate(SE=se_out2$SE, trait=phe)
  
  EffectMatrix[[t]] = eff_out2
}

EffectMatrix2 = bind_rows(EffectMatrix) 

###get the QTL effect for the significant markers identified for the full model
QTL_full_model = read.csv('QTLsWithFlankMarkers_full_model.csv')

QTL_Eff=NULL
for (tr in unique(QTL_full_model$lodcolumn)){
  tmp = QTL_full_model %>% filter(lodcolumn==tr)
  tmp2 = EffectMatrix2 %>% separate(trait,c('trait','Treatment'),sep = "_")%>% 
    filter(trait==tr)%>% filter(MARKER %in% tmp$MARKER) %>% left_join(QTL_full_model %>% dplyr::select(MARKER, QTL))
  QTL_Eff = rbind(QTL_Eff,tmp2)
}

QTL_Eff$Treatment = factor(QTL_Eff$Treatment, levels = c('subamb','ambient','elevated'))

pdf('../figures/QTL_Effect_Plots_ByTrait_ByMarker.pdf',onefile = T, width = 14, height = 7)
###plot the QTL effect
for(tr in unique(QTL_Eff$trait)){
    tmp = QTL_Eff %>% filter(trait==tr)
     print(
       ggplot(tmp, aes(Treatment, EFF)) + geom_point(size=2) + facet_grid(CROSS~QTL) + 
            ggtitle(tr)+
            geom_hline(yintercept = 0,linetype=2)+ xlab("")+theme_bw()+
            theme(axis.text = element_text(face = "bold", size = 15))+
            theme(strip.text = element_text(face="bold", size=15))+
            geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.1) + 
            ylab('QTL Effect') +xlab("") +coord_flip()
     )
}
dev.off()



