library(Seurat) # v5.0.3
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA) # v0.3.3
library(ggplot2)
library(ggsignif)
library(exactRankTests)

theme_set(theme_cowplot())
set.seed(1234)

packageVersion("Seurat")
packageVersion("hdWGCNA")

# Load a object
lung <- readRDS("lung_hdwgcna.rds")
lung

# Load subtype information
traitDATA = read.table("id_subtype_NE_nonNE.txt", check.names=FALSE, header=T) # please prepare a list of sample names and subtypes
NE_nonNE <- traitDATA$NE_nonNE
names(NE_nonNE) <- traitDATA$ID
OK <- NE_nonNE[lung$original]
names(OK) <- names(lung$original)
lung$NE_nonNE <- OK

lung$NE_nonNE <- factor(lung$NE_nonNE, levels = c("NE_ASCL1", "nonNE_ASCL2", "nonNE_YAP1", "Mixed", "DP"))

# DP vs Others
lung$DP <- as.vector(lung$NE_nonNE)
lung$DP[which(lung$NE_nonNE != "DP")] <- "Others"
DP <- lung$DP

group1 <- lung@meta.data %>% subset(DP == "DP") %>% rownames
length(group1)
group2 <- lung@meta.data %>% subset(DP == "Others") %>% rownames
length(group2)

# Differential ME analysis between DP and Others
DMEs <- FindDMEs(lung, barcodes1 = group1, barcodes2 = group2, test.use='wilcox', fc.slot = "counts")
write.table(DMEs, "DME_DP_vs_Others.txt",append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)


# Extraction of MEs for violin plots
MEs <- GetMEs(lung)
me <- colnames(MEs)
newMEs <- cbind(MEs, DP)

# Violin plots of MEs
pdf("violin_DP_vs_Others.pdf")
for(i in 1:ncol(MEs)) {
      p <- ggplot(newMEs, aes(x = DP, y = !!as.name(me[i]), fill = DP))
      p <- p + geom_violin()
      p <- p + geom_signif(comparisons = list(c("DP", "Others")), test = "wilcox.exact")
      p <- p + scale_fill_manual(values = c("purple", "grey"))
      print(p)
      P <- wilcox.exact(newMEs[,me[i]][which(newMEs$DP == "DP")], newMEs[,me[i]][which(newMEs$DP == "Others")])
      print(paste0(me[i],":p-value ", P$p.value))
      med <- median(newMEs[,me[i]][which(newMEs$DP == "DP")])
      print(paste0(me[i],":DP_median: ", med))
      med <- median(newMEs[,me[i]][which(newMEs$DP == "Others")])
      print(paste0(me[i],":Others_median: ", med))
}
dev.off()

# Violin plots of MEs (adjusted netgative values to zero)
# same as FindDMEs
MEs[MEs < 0] <- 0
me <- colnames(MEs)
newMEs <- cbind(MEs, DP)

pdf("violin_DP_vs_Others_round.pdf")
for(i in 1:ncol(MEs)) {
      p <- ggplot(newMEs, aes(x = DP, y = !!as.name(me[i]), fill = DP))
      p <- p + geom_violin()
      p <- p + geom_signif(comparisons = list(c("DP", "Others")), test = "wilcox.exact")
      p <- p + scale_fill_manual(values = c("purple", "grey"))
      print(p)
      P <- wilcox.exact(newMEs[,me[i]][which(newMEs$DP == "DP")], newMEs[,me[i]][which(newMEs$DP == "Others")])
      print(paste0(me[i],":p-value ", P$p.value))
      med <- log(mean(newMEs[,me[i]][which(newMEs$DP == "DP")]), 2) - log(mean(newMEs[,me[i]][which(newMEs$DP == "Others")]), 2)
      print(paste0(me[i],":avg.Log2FC: ", med))
}

dev.off()

rm(list = ls())
gc();gc()
