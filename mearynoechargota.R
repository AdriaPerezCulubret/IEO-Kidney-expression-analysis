library(SummarizedExperiment)
library(edgeR)

se <- readRDS(file.path("Documentos/MASTER/IEO/KIDNEY/seKIRC.rds"))
se

dim(colData(se))
clvar <- colnames(colData(se))
mcols(colData(se), use.names=TRUE)

rowRanges(se)

dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM


# MA PLOT NO NORM
dge$samples$group <- se$type
table(dge$samples$group)
plotSmear(dge, lowess=TRUE)


# MA PLOT NORM
dgenorm <- calcNormFactors(dge)
dgenorm$samples$group <- se$type
table(dgenorm$samples$group)
plotSmear(dgenorm, lowess=TRUE)

