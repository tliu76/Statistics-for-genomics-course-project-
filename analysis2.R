library(SummarizedExperiment)
library(matrixStats)


load("F:/2022 fourth term/statistic for genomics/rub_taybi_rnaseq_counts_matrix.rda")
se <- SummarizedExperiment(assays = SimpleList(exprs = features_by_samples_mat))
snames <- colnames(features_by_samples_mat)
snames.sp <- strsplit(snames, "_")

colData(se)$mouse <- paste0("m", sapply(snames.sp, "[[", 1))
colData(se)$celltype <- sapply(snames.sp, "[[", 2)
colData(se)$genotype <- sapply(snames.sp, "[[", 3)

b.se <- se[, se$celltype == "B"]
wt.se <- se[, se$genotype == "WT"]

b.se <- b.se[rowSums(assay(b.se)) > 12,]

libs <- colSums(assay(b.se)) / (30*10^6)

b.se3 <- b.se2 <- b.se
assay(b.se2) <- log2( assay(b.se) + 1)
assay(b.se3) <- log2( t(t(assay(b.se)) / libs) + 2)

pca1 <- prcomp(t(assay(b.se)))
summary(pca1)
par(mar=c(1,1,1,1))
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
dds <- DESeqDataSet(b.se, design = ~ genotype)
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds <- DESeq(dds)
res <- results(dds)

plotMA(res,)

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("genotype"))
BiocManager::install("sva")

library(sva)
dat <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat <- dat[idx,]

mod <- model.matrix( ~ genotype, colData(dds))
mod0 <- model.matrix( ~ 1, colData(dds))
svs <- svaseq(dat, mod, mod0)

## return to SVA
ddsva <- dds
ddsva$SV1 <- svs$sv[,1]
ddsva$SV2 <- svs$sv[,2]
design(ddsva) <- ~ genotype + SV1 + SV2
ddsva <- DESeq(ddsva)
resva <- results(ddsva, contrast = c("genotype", "CBP", "WT"))


## edgeR
library(edgeR)
dge <- DGEList(counts = assay(b.se), ref= "WT" ))
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
