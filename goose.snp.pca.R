library(ggplot2)
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)
library(gridExtra)
library(mypackage)

#liberal unlinked
goose.vcf <- read.vcfR("~/Downloads/biallelic.maf.unlinked.goose.recode.vcf")
#write out SNAPP formatted nexus file for building snapp tree
mypackage::vcf2SNAPP(goose.vcf, file = "~/Downloads/goose.nex")

#convert to genlight
goose.gen <- vcfR2genlight(goose.vcf, n.cores=1)
goose.gen@ind.names

#makepca1
pca <- glPca(goose.gen, nf=10)
pca.scores<-as.data.frame(pca$scores)
#ggplot color by species
ggplot(pca.scores, aes(x=PC1, y=PC2, col=as.factor(substr(goose.gen$ind.names, 1,9)))) +
  geom_point(cex = 2)+
  ggtitle("goose")+theme(legend.title = element_blank()) 

goose.gen@pop<-as.factor(substr(goose.gen$ind.names, 1,3))
neisd <- stamppNeisD(goose.gen, pop = FALSE)
plot(nj(neisd), type = "unrooted", cex = .5)

colnames(neisd) <- rownames(neisd)
heatmap.2(neisd, trace="none", cexRow=.8, cexCol=.8)




