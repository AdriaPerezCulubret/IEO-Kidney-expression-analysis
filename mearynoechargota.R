library(SummarizedExperiment)

se <- readRDS(file.path("Documentos/MASTER/IEO/KIDNEY/seKIRC.rds"))
se

dim(colData(se))

# colnames(colData(se))
