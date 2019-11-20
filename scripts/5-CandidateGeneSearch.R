setwd("~/Desktop/BranchyData/data/")
library(tidyverse)
library(ape)
options(warn = -1)

rm(list=ls())

### Read QTL info
QTL_all = read.csv('QTLsWithFlankMarkers_full_model.csv')

QTL = QTL_all %>% dplyr::mutate_at(vars(flank_lo, flank_hi, lodcolumn), list(as.character)) %>% 
  separate(flank_lo, c('chr','pos_lo')) %>% separate(flank_hi, c('chr_hi','pos_hi')) %>% 
  dplyr::mutate_at(vars(pos_lo, pos_hi), list(as.numeric))%>%
  dplyr:: select(chr,QTL, pos_lo, pos_hi, lodcolumn) %>% mutate(chr=str_replace(chr,"D",'d'))

gff3 = read.gff('Bdistachyon_314_v3.1.gene_exons.gff3.gz')
genes = gff3 %>% filter(type=='gene') %>% dplyr::rename(ID=attributes) %>%  separate(ID, c('A','locusName','C'), sep="=|\\;|\\.") %>%
 dplyr::select(seqid, start, end, locusName) 

annot = read.delim2('Bdistachyon_314_v3.1.annotation_info.txt',sep = "\t",header = T)
annot2 = annot %>% dplyr::select(locusName, GO, "Best.hit.arabi.name","arabi.symbol","arabi.defline","Best.hit.rice.name" , "rice.defline" )%>%
  separate(Best.hit.arabi.name,c('Best.hit.arabi.name','ver'),sep = '\\.') %>% separate(Best.hit.rice.name, c('Best.hit.rice.name','ver2'),sep = '\\.') %>%
  dplyr::select(-ver, -ver2) %>% distinct(locusName,.keep_all = T)

###write the annotation file for enrichment analysis if needed
geneID2GO = annot2 %>% dplyr::select(locusName, GO) %>% dplyr::filter(GO != "") 
write.table(geneID2GO, file ='annotation_for_enrichment.txt',sep=" ", quote = F, row.names = F, col.names = F )
####

genes_annot = genes %>% left_join(annot2) 

#######get the candidate genes for each QTL intervel for the traits
Counts = NULL
CandidateGeneAll = NULL
for (i in 1:nrow(QTL)){ 
  gene_list = genes_annot %>% dplyr::filter(seqid==QTL$chr[i]) %>% 
    filter(start >= as.numeric(QTL$pos_lo[i]) & end <= as.numeric(QTL$pos_hi[i]))
  
  QTL_sub=bind_rows(replicate(nrow(gene_list), QTL[i,], simplify = FALSE))
  gene_list = bind_cols(QTL_sub, gene_list)
  
  counts = data.frame(QTL=QTL$QTL[i], Trait=QTL$lodcolumn[i], BdGenes=length(gene_list$locusName), ArabGenes=length(unique(gene_list$Best.hit.arabi.name))-1,
                      RiceGenes=length(unique(gene_list$Best.hit.rice.name))-1)
  Counts = rbind(Counts, counts)
  write.csv(gene_list,paste0('CandidateGeneLists_',QTL$QTL[i],'_',QTL$lodcolumn[i], '.csv'), row.names = F)
  CandidateGeneAll = rbind(CandidateGeneAll,gene_list)
}
write.csv(Counts,'NumberOfCandidateGenesSummary.csv',row.names = F)
write.csv(CandidateGeneAll,'CandidateGene_ALL.csv',row.names = F)

###END