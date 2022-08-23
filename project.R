library(SummarizedExperiment)
library(matrixStats)
load("C:/Users/nizha/Downloads/hwkdata.rda")
se <- SummarizedExperiment(assays = SimpleList(exprs = features_by_samples_mat))
snames <- colnames(features_by_samples_mat)
snames.sp <- strsplit(snames, "_")

colData(se)$mouse <- paste0("m", sapply(snames.sp, "[[", 1))
colData(se)$celltype <- sapply(snames.sp, "[[", 2)
colData(se)$genotype <- sapply(snames.sp, "[[", 3)

b.se <- se[, se$celltype == "B"]
wt.se <- se[, se$genotype == "WT"]

wt.se <- wt.se[rowMeans(assay(wt.se),na.rm = T) > 1, ]

libs <- colSums(assay(wt.se)) / (30*10^6)

wt.se3 <- wt.se2 <- wt.se
assay(wt.se2) <- log2( assay(wt.se) + 1)
assay(wt.se3) <- log2( t(t(assay(wt.se)) / libs) + 2)

pca1 <- prcomp(t(assay(wt.se)))
summary(pca1)
plot(pca1$x[,1], pca1$x[,2], col = as.factor(wt.se$celltype), pch = 16,
     main = "counts")

pca2 <- prcomp(t(assay(wt.se2)))
summary(pca2)
plot(pca2$x[,1], pca2$x[,2], col = as.factor(wt.se2$celltype), pch = 16,
     main = "log2(counts + 1)")

pca3 <- prcomp(t(assay(wt.se3)), scale.=TRUE)
summary(pca3)
plot(pca3$x[,1], pca3$x[,2], col = as.factor(wt.se3$celltype), pch = 16,
     main = "log2(counts / depth + 1)")


libs = colSums(assay(wt.se))

## MA PLOT
a = rowMeans(log2(    t(t(assay(wt.se))/libs)[,1:2] +1))
m = log2(t(t(assay(wt.se))/libs)+1)[,1] - log2(t(t(assay(wt.se))/libs)+1)[,2]
library(scales)
plot(a,m, pch = 16, cex = 0.5, col = alpha("black", 0.5))

## DEseq2
library(DESeq2)
assay(wt.se) <- round(assay(wt.se))
dds <- DESeqDataSet(wt.se, design = ~ celltype)
dds$celltype <- relevel(dds$celltype, ref = "B")
dds <- DESeq(dds)
res <- results(dds)
plotMA(res, main = "MA plot comparing differential expressions in T cells vs. B cells (ref)", ylim = c(-10,10))
res[res$padj < 0.01 & !is.na(res$padj),]
res[res$baseMean<0.1,]
sum(res$padj < 0.1, na.rm = T)
hist(res$pvalue, main = "Histogram of p-values", xlab = "p-values")
wt.se

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("celltype"))

res_sig <- subset(res, padj<.1)
res_sig_sorted = res_sig[order(res_sig$padj), ]
head(res_sig_sorted,n = 10)

library(tidyverse)
library(sva)
dat <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat <- dat[idx,]

mod <- model.matrix( ~ celltype, colData(dds))
mod0 <- model.matrix( ~ 1, colData(dds))
svs <- svaseq(dat, mod, mod0)


## return to SVA
ddsva <- dds
ddsva$SV1 <- svs$sv[,1]
ddsva$SV2 <- svs$sv[,2]
ddsva$SV3 <- svs$sv[,3]
design(ddsva) <- ~ celltype + SV1 + SV2 + SV3
ddsva <- DESeq(ddsva)
resva <- results(ddsva, contrast = c("celltype", "T", "B"))
resva
sum(resva$padj < 0.1, na.rm = T)
resva.sig = resva[resva$padj<0.1,]

############plotting
plotMA(resva, main = "MA plot of surrogate variable analysis comparing differential expressions in T cells vs. B cells (ref)", ylim = c(-10,10))
plot1 = as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(plot1) = c("FDR Threshold","Number of Significant Genes")
plot1[,1] = 10^(c(-1,-2,-5,-10,-20))
plot1[,2] = c(sum(res_sig$padj<0.1),sum(res_sig$padj<0.01),sum(res_sig$padj<1*10^-5),sum(res_sig$padj<1*10^-10),sum(res_sig$padj<1*10^-20))


plot2 = as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(plot2) = c("FDR Threshold","Number of Significant Genes")
plot2[,1] = plot1[,1]
plot2[,2] = c(sum(resva.sig$padj<0.1),sum(resva.sig$padj<0.01),sum(resva.sig$padj<1*10^-5),sum(resva.sig$padj<1*10^-10),sum(resva.sig$padj<1*10^-20))
plot1$method = "unadjusted for SV"
plot2$method = "adjusted for SVs"
wt.se_result = rbind(plot1,plot2)

wt.se_result %>%
  ggplot(aes(x = -log10(`FDR Threshold`), y = `Number of Significant Genes`, color = method))+
  geom_point()+
  geom_text(aes(label = `Number of Significant Genes`), nudge_x = T, color = "blue")+
  geom_line()+
  theme_bw()

## edgeR
library(edgeR)
dge <- DGEList(counts = assay(wt.se), group = relevel(factor(wt.se$celltype),ref = "B"))
keep <- filterByExpr(dge)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log=TRUE)
pca <- prcomp(t(logcpm), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col = as.factor(wt.se$celltype), pch = 16)

design <- model.matrix(~ group, data = dge$samples)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt)

library(sva)
dat <- cpm(dge)
idx  <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix( ~ group, data = dge$samples)
mod0 <- model.matrix( ~ 1, data = dge$samples)
svs2 <- svaseq(dat, mod, mod0)


dge2 <- dge
dge2$samples$SV1 <- svs2$sv[,1]
dge2$samples$SV2 <- svs2$sv[,2]
dge2$samples$SV3 <- svs2$sv[,3]
design2 <- model.matrix(~ group + SV1 + SV2 + SV3, data = dge2$samples)
dge2 <- estimateDisp(dge2, design2)
fit2 <- glmFit(dge2, design2)
lrt2 <- glmLRT(fit2, coef = 2)
tt2 <- topTags(lrt2)

#####################b cell
b.se <- se[, se$celltype == "B"]

b.se <- b.se[rowSums(assay(b.se)) > 12,]

libs <- colSums(assay(b.se)) / (30*10^6)

b.se3 <- b.se2 <- b.se
assay(b.se2) <- log2( assay(b.se) + 1)
assay(b.se3) <- log2( t(t(assay(b.se)) / libs) + 2)

pca1 <- prcomp(t(assay(b.se)))
summary(pca1)
plot(pca1$x[,1], pca1$x[,2], col = as.factor(b.se$genotype), pch = 16,
     main = "counts")

pca2 <- prcomp(t(assay(b.se2)))
summary(pca2)
plot(pca2$x[,1], pca2$x[,2], col = as.factor(b.se2$genotype), pch = 16,
     main = "log2(counts + 1)")

pca3 <- prcomp(t(assay(b.se3)), scale.=TRUE)
summary(pca3)
plot(pca3$x[,1], pca3$x[,2], col = as.factor(b.se3$genotype), pch = 16,
     main = "log2(counts / depth + 1)")




## DEseq2
library(DESeq2)
assay(b.se) <- round(assay(b.se))
dds.b <- DESeqDataSet(b.se, design = ~ genotype)
dds.b$genotype <- relevel(dds.b$genotype, ref = "WT")
dds.b <- DESeq(dds.b)
res.b <- results(dds.b)
plotMA(res.b, ylim = c(-5,5))
sum(res.b$padj < 0.01, na.rm = T)
vsd.b <- vst(dds.b, blind = TRUE)
plotPCA(vsd.b, intgroup = c("genotype"))
res.b_sig = subset(res.b, res.b$padj<0.1)

library(sva)
dat.b <- counts(dds.b, normalized = TRUE)
idx.b  <- rowMeans(dat.b) > 1
dat.b <- dat.b[idx.b,]

mod.b <- model.matrix( ~ genotype, colData(dds.b))
mod0.b <- model.matrix( ~ 1, colData(dds.b))
svs.b <- svaseq(dat.b, mod.b, mod0.b)

## return to SVA
ddsva.b <- dds.b
ddsva.b$SV1 <- svs.b$sv[,1]
ddsva.b$SV2 <- svs.b$sv[,2]
design(ddsva.b) <- ~ genotype + SV1 + SV2
ddsva.b <- DESeq(ddsva.b)
resva.b <- results(ddsva.b, contrast = c("genotype", "CBP", "WT"))
resva.b_sig = subset(resva.b, resva.b$padj<0.1)

############plotting
plot3 = as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(plot3) = c("FDR Threshold","Number of Significant Genes")
plot3[,1] = 10^(c(-1,-2,-5,-10,-20))
plot3[,2] = c(sum(res.b_sig$padj<0.1),sum(res.b_sig$padj<0.01),sum(res.b_sig$padj<1*10^-5),sum(res.b_sig$padj<1*10^-10),sum(res.b_sig$padj<1*10^-20))


plot4 = as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(plot4) = c("FDR Threshold","Number of Significant Genes")
plot4[,1] = plot3[,1]
plot4[,2] = c(sum(resva.b_sig$padj<0.1),sum(resva.b_sig$padj<0.01),sum(resva.b_sig$padj<1*10^-5),sum(resva.b_sig$padj<1*10^-10),sum(resva.b_sig$padj<1*10^-20))
plot3$method = "unadjusted for SV"
plot4$method = "adjusted for SVs"
b.se_result = rbind(plot3,plot4)

b.se_result %>%
  ggplot(aes(x = -log10(`FDR Threshold`), y = `Number of Significant Genes`, color = method))+
  geom_point()+
  geom_text(aes(label = `Number of Significant Genes`), nudge_x = T,nudge_y = T, color = "blue", check_overlap = T)+
  geom_line()+
  theme_bw()

## edgeR
library(edgeR)
dge <- DGEList(counts = assay(b.se),
               ref= "WT"))
keep <- filterByExpr(dge)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log=TRUE)
pca <- prcomp(t(logcpm), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], col = as.factor(b.se$genotype), pch = 16)

design <- model.matrix(~ group, data = dge$samples)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt)

library(sva)
dat <- cpm(dge)
idx  <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix( ~ group, data = dge$samples)
mod0 <- model.matrix( ~ 1, data = dge$samples)
svs2 <- svaseq(dat, mod, mod0)


dge2 <- dge
dge2$samples$SV1 <- svs2$sv[,1]
dge2$samples$SV2 <- svs2$sv[,2]
dge2$samples$SV3 <- svs2$sv[,3]
design2 <- model.matrix(~ group + SV1 + SV2 + SV3, data = dge2$samples)
dge2 <- estimateDisp(dge2, design2)
fit2 <- glmFit(dge2, design2)
lrt2 <- glmLRT(fit2, coef = 2)
tt2 <- topTags(lrt2)

####facet
wt.se_result$group = "T cells vs. B cells in Wildtype"
b.se_result$group = "mutant vs wildtype in B cells"
res_tot = rbind(wt.se_result,b.se_result)
res_tot %>%
  ggplot(aes(x = -log10(`FDR Threshold`), y = `Number of Significant Genes`, color = method))+
  geom_point()+
  geom_text(aes(label = `Number of Significant Genes`), nudge_x = T, color = "blue", check_overlap = T)+
  geom_line()+
  facet_grid(.~group, scales = "free_y")+
  theme_bw()