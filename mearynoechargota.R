library(SummarizedExperiment)
se <- readRDS(file.path("~/GitHub/KIDNEY/seKIRC.rds"))
se
dim(colData(se))
mcols(colData(se), use.names=TRUE)
rowRanges(se)
library(edgeR)
dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
dge
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[1:5, 1:5]
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
names(sampledepth) <- substr(colnames(se), 6, 12)
sort(sampledepth)

# GENE EXPRESSION LogCPM PLOT
library(geneplotter)
par(mfrow=c(1,2), mar=c(4,5,1,1))
# ALL SAMPLES PLOT
multidensity(as.list(as.data.frame(assays(se)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
# NORMAL AND TUMOR PLOTS
normalCPM <- se[,se$type == "normal"]
tumorCPM <- se[,se$type == "tumor"]
multidensity(as.list(as.data.frame(assays(normalCPM)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(tumorCPM)$logCPM)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# NORMALIZATION PLOTS
dgenorm <- calcNormFactors(dge)
assays(se)$logCPMnorm <- cpm(dgenorm, log=TRUE, prior.count=0.5)
normalCPMnorm <- se[,se$type == "normal"]
tumorCPMnorm <- se[,se$type == "tumor"]
multidensity(as.list(as.data.frame(assays(normalCPMnorm)$logCPMnorm)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(tumorCPMnorm)$logCPMnorm)), xlab = "log2 CPM", legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
