library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(geneplotter)

se <- readRDS(file.path("~/GitHub/KIDNEY/seKIRC.rds"))

clvar <- colnames(colData(se))
mcols(colData(se), use.names=TRUE)

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


# GENE EXPRESSION LogCPM PLOT

par(mfrow=c(1,2), mar=c(4,5,1,1))
# ALL SAMPLES PLOT
multidensity(as.list(as.data.frame(assays(se)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
# NORMAL AND TUMOR PLOTS
normalCPM <- se[,se$type == "normal"]
tumorCPM <- se[,se$type == "tumor"]
multidensity(as.list(as.data.frame(assays(normalCPM)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(tumorCPM)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# NORMALIZATION PLOTS
normalCPMnorm <- se[,se$type == "normal"]
tumorCPMnorm <- se[,se$type == "tumor"]
multidensity(as.list(as.data.frame(assays(normalCPMnorm)$logCPMnorm)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(tumorCPMnorm)$logCPMnorm)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
