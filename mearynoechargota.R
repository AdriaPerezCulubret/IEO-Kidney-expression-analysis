library(SummarizedExperiment)
library(edgeR)
library(ggplot2)

se <- readRDS(file.path("Documentos/MASTER/IEO/KIDNEY/seKIRC.rds"))
se

dim(colData(se))
clvar <- colnames(colData(se))
mcols(colData(se), use.names=TRUE)

rowRanges(se)

dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)


# NORMALIZATION
dgenorm <- calcNormFactors(dge)
assays(se)$logCPMnorm <- cpm(dgenorm, log=TRUE, prior.count=0.5)


# MA PLOT NO NORM
dge$samples$group <- se$type
table(dge$samples$group)
plotSmear(dge, lowess=TRUE)


# MA PLOT NORM
dgenorm$samples$group <- se$type
table(dgenorm$samples$group)
plotSmear(dgenorm, lowess=TRUE)


# GENE DIST

# NO NORM
genemean <- data.frame(Means=rowMeans(assays(se)$logCPM[,-1]))

ggplot(genemean) + 
  geom_bar(aes(x=Means), binwidth = 1.5, fill="white", color="black") +
  theme_bw() + xlab("\nlog2CPM") + ylab("Count\n") +
  geom_vline(xintercept=mean(genemean$Means), color="red")

# NORM

genemean.norm <- data.frame(Means=rowMeans(assays(se)$logCPMnorm[,-1]))

ggplot(genemean.norm) + 
  geom_bar(aes(x=Means), binwidth = 1.5, fill="white", color="black") +
  theme_bw() + xlab("\nlog2CPM") + ylab("Count\n") +
  geom_vline(xintercept=mean(genemean.norm$Means), color="red")

