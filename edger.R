# Fonte: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# carregar a biblioteca edgeR
library("edgeR")

# table with gene count
# First Column: Symbol
# The remaining columns are the gene count
TabelaCount <- "merge-table-STAR-gene-level-5-50x.csv"

# loading the table (the first column is called Symbol)
x <- read.delim(TabelaCount,row.names="symbol")

# groups: 1 and 2
# The first 4 samples group 1
# 5 afterwards is group 2
group <- factor(c(1,1,1,1,2,2,2,2,2))
y <- DGEList(counts=x,group=group)


# We filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# TMM normalization and display the normalization factors
y <- calcNormFactors(y)

png("plotMDS.png")
plotMDS(y)
dev.off()

design <- model.matrix(~group)
y <- estimateDisp(y,design)

png("plotBCV.png")
plotBCV(y)
dev.off()

#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

png("plotMD.png")
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

write.table(lrt$table,"lrt.DGElist.csv",quote=FALSE, sep="\t")
