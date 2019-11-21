setwd("~/Desktop/BranchyData/data/")
library(qtl)
#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtlTools)
library(tidyverse)
options(warn = -1)

rm(list = ls())
par(mfrow=c(1,1))
cross2 = read_cross2("RIL_Multiple.yaml")
map = insert_pseudomarkers(cross2$gmap, step=1)
pr_multiple = calc_genoprob(cross2, map, error_prob=0.002,cores = 2)
kinship_loco_multiple= calc_kinship(pr_multiple, "loco", cores=2)

covar = as.data.frame(cross2$covar)
covar$Treatment = as.numeric(covar$Treatment)

covar_factor = factor(covar$Treatment)
covar_matrix = model.matrix( ~ covar_factor)[ , -1]

covar = as.matrix(covar)
###full model: G + E + GXE, for all traits
out_full = scan1(pr_multiple, cross2$pheno, kinship_loco_multiple,addcovar = covar_matrix,intcovar = covar_matrix, cores = 2)
###permutation test with stratifier of covariates
set.seed(1234)
perm_full = scan1perm(pr_multiple, cross2$pheno, kinship_loco_multiple, addcovar = covar_matrix,
                       intcovar = covar_matrix,n_perm=1000,perm_strata = covar, core=3)

###reduced model: G + E
out_reduced = scan1(pr_multiple, cross2$pheno, kinship_loco_multiple,addcovar = covar_matrix, cores = 2)
###permutation test
set.seed(1234)
perm_reduced <- scan1perm(pr_multiple, cross2$pheno, kinship_loco_multiple, addcovar = covar_matrix,
                           n_perm=1000, perm_strata = covar, cores = 3)

save(out_full, perm_full, out_reduced, perm_reduced, file='scan1Results.RData')

###Get the threshold value and test GXE
gxe = out_full - out_reduced

threshold_full = summary(perm_full,0.05)
threshold_reduced = summary(perm_reduced,0.05)
threshold_gxe = threshold_full - threshold_reduced ## or use: threshold_gxe = summary(gxe,0.05)

###get the peak and flank markers for full model and gxe model
peaks_full = find_peaks(out_full, map, threshold=threshold_full, drop=1, peakdrop = 1.5) 
peaks_gxe = find_peaks(gxe, map, threshold=threshold_gxe, drop=1, peakdrop = 1.5) 

##merge the marker information into the significant QTLs
nnmap = data.frame(pos=unlist(map))
nnmap = nnmap %>% rownames_to_column('MARKER') %>% rowwise() %>% mutate(chr = str_sub(MARKER, 1,1), MARKER=str_sub(MARKER, 3,nchar(MARKER)))

peaks_full =peaks_full %>% left_join(nnmap)
peaks_gxe=peaks_gxe %>% left_join(nnmap)


flank_marker_full = peaks_full %>% mutate(flank_lo = find_marker(map, chr, ci_lo), 
                                          flank_hi=find_marker(map, chr, ci_hi)) %>%
                  mutate(QTL=paste0(chr, '@',round(pos, 1)))
flank_marker_gxe = peaks_gxe %>% mutate(flank_lo = find_marker(map, chr, ci_lo), 
                                        flank_hi=find_marker(map, chr, ci_hi)) %>%
                 mutate(QTL=paste0(chr, '@',round(pos, 1)))

write.csv(flank_marker_full,'QTLsWithFlankMarkers_full_model.csv') ## for later QTL effect plot
write.csv(flank_marker_gxe,'QTLsWithFlankMarkers_gxe_model.csv')

#####Plot the QTLs on genetic map
pdf('../figures/QTLs on Genetic Map.pdf', width = 9,height = 7, onefile = T)

cross = read.cross("csv",file="Brachypodium_RQTL data.csv")
cols  =  c('darkblue','blue',"green","red","violet")###add your colors depending on how many traits you have
with(flank_marker_full, segmentsOnMap(cross, phe = lodcolumn, chr = chr,
                                      l = ci_lo, h = ci_hi,
                                      peaklod = lod, peakcM = pos,
                                      showPeaks = TRUE,
                                      chrBuffer = c(0.15,0.1),
                                      tick.width=0.05, lwd=2,
                                      col=cols,
                                      leg.inset=0.05, legendCex=1,
                                      legendPosition="bottomright")
)
mtext('QTL Summary from the Full Model',side = 3)

with(flank_marker_gxe, segmentsOnMap(cross, phe = lodcolumn, chr = chr,
                                      l = ci_lo, h = ci_hi,
                                      peaklod = lod, peakcM = pos,
                                      showPeaks = TRUE,
                                      chrBuffer = c(0.15,0.1),
                                      tick.width=0.05, lwd=2,
                                      col=cols,
                                      leg.inset=0.05, legendCex=1,
                                      legendPosition="bottomright")
)
mtext('QTL Summary from the G x E Model',side = 3)
dev.off()

###plotting in diffrent ways
pdf('../figures/LOD Profiles for Different Traits with Different Models.pdf',width = 10,height = 7)
par(mfrow=c(5,1),oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,0.6,1) + 0.1)

ymx=max(maxlod(out_full))

for (i in 1:5){
plot(out_full, lodcolumn = i, map, col='red', ylim=c(0,ymx+0.2))
plot(out_reduced,lodcolumn = i, map, col='blue',add=T)
plot(gxe,lodcolumn = i, map, col='green',add=T)
abline(h =threshold_full[i], col='red', lty=2)
abline(h =threshold_reduced[i], col='blue', lty=2)
abline(h =threshold_gxe[i], col='green', lty=2)
legend("topleft",legend = colnames(out_full)[i],bty = "n")
}
mtext("LOD Score", side = 2, outer = T, line = 1.5)

legend(x=200,y=-0.5, legend=c('Full Model','Reduced Model','G x E'),cex = 2,
       lty = c(1,1,1),lwd = 2, col = c('red','blue','green'), horiz = T, bty="n", xpd=NA)

dev.off()

##########plot the traits together for full model and GXE model
pdf('../figures/LOD Profiles for All Traits with Full Models and GxE model.pdf',width = 10,height = 7 )
par(mfrow=c(2,1),oma = c(5,4,4,0) + 0.1,
    mar = c(0,0,0,1) + 0.1)
cols = c("darkblue","blue", "violetred","red","green" ) ###pick your color..

plot(out_full, map, lodcolumn = 1, col=cols[1], ylim= c(0,ymx+0.2), lty=1)
abline(h=threshold_full[1], col=cols[1], lty=2)

for (i in 2:(ncol(out_full))){
  plot(out_full, map, lodcolumn = i, col=cols[i], add=T)
  abline(h=threshold_full[i], col=cols[i], lty=2)
}
legend('topleft',"Full Model", bty="n")

plot(gxe, map, lodcolumn = 1, col=cols[1], ylim= c(0,ymx+0.2), lty=1)
abline(h=threshold_gxe[1], col=cols[1], lty=2)

for (i in 2:(ncol(gxe))){
  plot(gxe, map, lodcolumn = i, col=cols[i], add=T)
  abline(h=threshold_gxe[i], col=cols[i], lty=2)
}
legend('topleft',"G x E Model", bty = "n")

mtext("LOD Score", side = 2, outer = T, line = 1.5)
legend(x=200,y=-0.5, legend=colnames(out_full),cex = 1.5,
       lwd=2, col = cols, bty="n", xpd=NA,  ncol = 3)
dev.off()
####END

