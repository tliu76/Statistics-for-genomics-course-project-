library(SummarizedExperiment)
setwd("F:/2022 fourth term/statistic for genomics")
load("F:/2022 fourth term/statistic for genomics/rub_taybi_rnaseq_counts_matrix.rda")
se = SummarizedExperiment(counts = features_by_samples_mat)
se = SummarizedExperiment(assays  = SimpleList(counts = features_by_samples_mat))
cnames.sp = strsplit(colnames(se), "_")
colData(se)$celltype = sapply(cnames.sp, function(xx) xx[2])
colData(se)$genotype = sapply(cnames.sp, function(xx) xx[3])
colData(se)$mouse = paste0("m", sapply(cnames.sp, function(xx) xx[1]))
colData(se)
se.b = se[, se$celltype == "B"]
libs = colSums(assay(se))

## MA PLOT
a = rowMeans(log2( t(t(assay(se.b))/libs)[,1:2] +1))
m = log2(t(t(assay(se.b))/libs)+1)[,1] - log2(t(t(assay(se.b))/libs)+1)[,2]
library(scales)
plot(a,m, pch = 16, cex = 0.5, col = alpha("black", 0.05))

## PCA 
se.b  = se.b[rowSums(assay(se.b)) > 12, ]
pca = prcomp(t(log2(t(t(assay(se.b))/ libs+1))), scale. = TRUE)
summary(pca)
plot(pca$x[,1], pca$x[,2], col = as.factor(se.b$genotype), pch = 16)
