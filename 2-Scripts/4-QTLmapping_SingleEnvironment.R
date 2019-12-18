setwd("~/Desktop/BranchyData/1-ProcessedData/")
library(qtl)
library(sommer)
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(tidyverse)
options(warn = -1)

rm(list=ls())

###########Heritability and variance, genetic correlation using Sommer package
##1. format the files and get the Additive kinship matrix
geno_raw = read.csv('RIL_Single_gen.csv')
geno = geno_raw[,-1]
geno  = geno %>% mutate_if(.,is.factor, as.character) %>% 
  mutate_all(funs(str_replace(.,"A", "1"))) %>% mutate_all(funs(str_replace(.,"B", "-1"))) %>%
  mutate_if(., is.character, as.numeric)
rownames(geno) = geno_raw$id
geno = data.matrix(geno)

A = A.mat(geno)

###read in the phenotype data
phenos = read.csv('Phenotypes_LinedUpWithGenotypes.csv')

pheid = colnames(phenos)[-c(1:2)]
heritability = c()
varComp = c()

for (s in unique(phenos$Treatment)){
  df1 = phenos %>% dplyr::filter(Treatment==s) 
  for (phe in pheid){
    print(c(s,phe))
    if (sum(df1[,phe],na.rm = T)==0) {h2 = data.frame(Treatment=s,TRAIT=phe,h2=NA,h2_SE=NA)
    Vave = cbind(data.frame(Treatment=rep(s,2),TRAIT=rep(phe,2),VARIANCE=c('Va','Ve'), 
                            VarComp=rep(NA,2),VarCompSE=rep(NA,2),Zratio=rep(NA,2)))}else{
                              form =as.formula(paste(phe," ~ 1"))
                              A_phe = mmer(form, random=~vs(genotype, Gu=A), rcov=~vs(units),data=df1, verbose = F)
                              VaVe = cbind(data.frame(Treatment=rep(s,2), TRAIT=rep(phe,2),VARIANCE=c('Va','Ve')),data.frame(summary(A_phe)$varcomp))
                              h2 = data.frame(c(s, phe,round(pin(A_phe, h2 ~ V1 / ( V1 + V2 )),2)))
                              colnames(h2) = c('Treatment','TRAIT','h2','h2_SE')}
    heritability = rbind(heritability,h2)
    varComp = rbind(varComp,VaVe)}
}
write.csv(heritability, "HeritabilityForEachTraitAtEachTreatment_usingSommer.csv",row.names = F)


###plot heritability
pdf("../3-Figures/3-HeritabilityVariance.pdf", onefile = T, width = 12, height = 8)
ggplot(heritability, aes(Treatment,h2)) + geom_errorbar(aes(ymin=0.005, ymax=h2+h2_SE), width=0.2)+ 
  geom_col() + facet_wrap(TRAIT~., nrow = 3)+
  theme_bw()+ ylab('Heritability')+xlab('')+
  theme(text = element_text(face = "bold", size = 12))+
  theme(axis.text = element_text(face = "bold", size = 12),axis.title = element_text(face = "bold", size = 14))+
  theme(strip.text = element_text(face="bold", size=12))


###plot variance components
vc = varComp
vc_va = subset(vc, VARIANCE=='Va')
vc_ve = subset(vc, VARIANCE=='Ve')
vc_ve = vc_ve %>% dplyr::mutate(VarComp=-VarComp)

ggplot(vc_va,aes(Treatment,VarComp,fill=VARIANCE))+geom_col() + facet_wrap(TRAIT~., nrow = 3, scales = "free")+
  geom_hline(yintercept=0, color = "black",lwd=1.2)+ylab('')+xlab('')+
  geom_col(data=vc_ve,aes(Treatment,VarComp,fill=VARIANCE))+theme_bw()

###plot genetic correlation
##XXXXXX
dev.off()

###genetic correlation for each trait among treatments
#pdf("../3-Figures/3-GeneticCorrelation.pdf", height = 8, width = 10)
Rg <- c()
for (i in 3:ncol(phenos)){
  tmp = phenos[,c(1,2,i)] 
  colnames(tmp)[3] = "Trait"
  tmp = tmp %>% spread(Treatment,Trait)
 # print(colnames(phenos)[i])
  rg1 = mmer(cbind(subamb,ambient,elevated)~1, random=~vs(genotype, Gu=A), rcov=~units,data=tmp, tolparinv = 1e-02) 
  x1 = as.data.frame(cov2cor(rg1$sigma$`u:genotype`)) %>% 
    rownames_to_column('Treatment') %>% dplyr::mutate(Trait=colnames(phenos)[i])
  Rg = rbind(Rg,x1)
}

write.csv(Rg,'GeneticCorrelation_each_site.csv',row.names=F)



######running Rqtl2 for each trait at each treatment
###using Rqtl2 for RIL population
cross2 = read_cross2("RIL_Single.yaml")
map = insert_pseudomarkers(cross2$gmap, step=1)
pr_single = calc_genoprob(cross2, map, error_prob=0.002,cores = 2)

###get the h2 and va/ve using kinship matrix
####
kinship_all_single = calc_kinship(pr_single, core=2)
pheno = cross2$pheno
h2 = est_herit(pheno,kinship_all_single)
h2 = as.data.frame(h2)
write.csv(h2, "HeritabilityForEachTraitAtEachTreatment_usingRqtl2.csv",row.names = T)

###plot the h2 from Rqtl2

h2_all = h2  %>%  rownames_to_column("id")%>% filter(str_detect(id, "ambient|subamb|elevated"))%>%
  separate(id, c("trait","A1","A2","Treatment"), fill="left") %>%
  unite("TRAIT", trait, A1, A2) %>% mutate_at(vars(h2), round,2)%>%
  mutate_all(funs(str_replace(.,"NA_", ""))) %>% mutate_all(funs(str_replace(.,"NA_", ""))) %>% 
  filter(complete.cases(.))%>% mutate_at(vars(h2), as.numeric)

h2_all$TRAIT = factor(h2_all$TRAIT, levels = c("biomass","biomass_log","leaf_mass","leaf_mass_log",
                                               "seed_mass","seed_mass_log","stem_mass","stem_mass_log",
                                               "ht_flag","ht_flower","AmountC","AmountN","d13C","d15N","CN"))
h2_all$Treatment = factor(h2_all$Treatment, levels= c("subamb","ambient","elevated"))

pdf("../3-Figures/3-Heritability_Rqtl2.pdf", onefile = T, width = 12, height = 8)
ggplot(h2_all, aes(Treatment,h2)) + 
  geom_col() + facet_wrap(TRAIT~., nrow = 3)+
  theme_bw()+ ylab('Heritability')+xlab('')+
  theme(text = element_text(face = "bold", size = 12))+
  theme(axis.text = element_text(face = "bold", size = 12),axis.title = element_text(face = "bold", size = 14))+
  theme(strip.text = element_text(face="bold", size=12))

dev.off()




kinship_loco_single = calc_kinship(pr_single, "loco", cores=2)

###run scan1 
out = scan1(pr_single, pheno = pheno[,62:69], kinship_loco_single)
out_perm = scan1perm(pr_single,pheno = pheno[,62:69], kinship_loco_single, n_perm = 1000, cores = 3)
thresholds = summary(out_perm,0.05)
##
pdf("../3-Figures/4-LOD profile for averages across CO2 and differnces betewen CO2.pdf",onefile = T, width = 10,height = 7)

cols  = c("brown","red", RColorBrewer::brewer.pal(n=11,name = "BrBG"),"darkgrey","black")
#cols = c("darkblue","blue", "violetred","red","green" ) ###pick your color..
plot(out,lodcolumn = 70, map, col=cols[1], ylim=c(0,4))
abline(h=thresholds[70], col=cols[1], lty=2)
for (i in 71:84){
  plot(out, map, lodcolumn = i, col=cols[i-69],add=T)
  abline(h=thresholds[i], col=cols[i-69], lty=2)
}
legend('topleft',legend = colnames(out)[70:84],col=cols, lty = 1, bty="n", ncol = 5)
mtext("QTL mapping of traits across CO2 levels", outer = F, side=3)


##the differnces between CO2 levels
high_level = grep("high[^0-9]",colnames(out))
plot(out,lodcolumn = high_level[1], map, col=cols[1], ylim=c(0, 4))
abline(h=thresholds[high_level[1]], col=cols[1], lty=2)

for (i in high_level[-1]){
  plot(out, map, lodcolumn = i, col=cols[(i-42)/2],add=T)
  abline(h=thresholds[i], col=cols[(i-42)/2], lty=2)
}
legend('topleft',legend = colnames(out)[high_level],col=cols, lty = 1, bty="n", ncol = 4)
mtext("QTL mapping of traits between CO2 elevated and ambient", outer = F, side=3)

low_level = grep("_low",colnames(out))
plot(out,lodcolumn = low_level[1], map, col=cols[1], ylim=c(0, 4))
abline(h=thresholds[low_level[1]], col=cols[1], lty=2)
for (i in low_level[-1]){
  plot(out, map, lodcolumn = i, col=cols[(i-43)/2],add=T)
  abline(h=thresholds[i], col=cols[(i-43)/2], lty=2)
}
legend('topleft',legend = colnames(out)[low_level],col=cols, lty = 1, bty="n", ncol=4)
mtext("QTL mapping of traits between CO2 subambient and ambient", outer = F, side=3)


high2_level = grep("high2",colnames(out))
plot(out,lodcolumn = high2_level[1], map, col=cols[1], ylim=c(0, 4))
abline(h=thresholds[high2_level[1]], col=cols[1], lty=2)

for (i in high2_level[-1]){
  plot(out, map, lodcolumn = i, col=cols[i],add=T)
  abline(h=thresholds[i], col=cols[i], lty=2)
}
legend('topleft',legend = colnames(out)[high2_level],col=cols, lty = 1, bty="n", ncol = 4)
mtext("QTL mapping of traits between CO2 elevated and subamb", outer = F, side=3)


dev.off()

pheid = colnames(cross2$pheno)[1:45]
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
EffectMatrix3 = EffectMatrix2 %>%  separate(trait, c("trait","A1","A2","Treatment"), fill="left") %>%
  unite("trait", trait, A1, A2)%>% mutate_all(funs(str_replace(.,"NA_", ""))) %>% mutate_all(funs(str_replace(.,"NA_", "")))
###get the QTL effect for the significant markers identified for the full model
QTL_full_model = read.csv('QTLsWithFlankMarkers_full_model.csv')

QTL_Eff=NULL
for (tr in unique(QTL_full_model$lodcolumn)){
  tmp = QTL_full_model %>% filter(lodcolumn==tr)
  tmp2 = EffectMatrix3%>% filter(trait==tr)%>% filter(MARKER %in% tmp$MARKER) %>% left_join(QTL_full_model %>% dplyr::select(MARKER, QTL))
  QTL_Eff = rbind(QTL_Eff,tmp2)
}

QTL_Eff$Treatment = factor(QTL_Eff$Treatment, levels = c('subamb','ambient','elevated'))
QTL_Eff$EFF = as.numeric(QTL_Eff$EFF)
QTL_Eff$SE = as.numeric(QTL_Eff$SE)

pdf('../3-Figures/4-QTL_Effect_Plots_ByTrait_ByMarker.pdf',onefile = T, width = 14, height = 7)
###plot the QTL effect
for(tr in unique(QTL_Eff$trait)){
    tmp = QTL_Eff %>% filter(trait==tr)
     print(
       ggplot(tmp, aes(Treatment, EFF, fill=CROSS)) + geom_col(size=2) + facet_grid(.~QTL) + 
            ggtitle(tr)+
            geom_hline(yintercept = 0,linetype=2)+ xlab("")+theme_bw()+
            theme(axis.text = element_text(face = "bold", size = 15))+
            theme(strip.text = element_text(face="bold", size=15))+
            geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), width=0.1, position=position_dodge(0.05) ) + 
            ylab('QTL Effect') +xlab("")
     )
}
dev.off()



