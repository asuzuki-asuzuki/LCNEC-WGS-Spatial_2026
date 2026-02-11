library(gplots)

set.seed(1234)

# Load a expression matrix of module 17
DATA <- read.table("exp_module17.txt", header=1, row.names = 1)
DATA2 <- log2(DATA + 1)

# Clustering
d1 <- dist(DATA2, method="euclid")
d2 <- dist(t(DATA2), method="euclid")
c1 <- hclust(d1, method = "ward.D2")
c2 <- hclust(d2, method = "ward.D2")

# Extract 4 clusters
clusters <- cutree(c2, k = 4)
col <- c("blue", "red", "orange", "lightblue")[clusters]

# Z-score
DATA3 <- t(scale(t(DATA2)))

# Heatmap
pdf("heatmap_Module17.pdf", width = 9, height = 7)
heatmap.2(as.matrix(DATA3), Colv=as.dendrogram(c2), Rowv=as.dendrogram(c1), dendrogram = "both", scale = "none", trace="none", density.info="none",
  col=colorpanel(44, low="magenta", mid="black", high="green"), keysize = 1, symkey=F, margins=c(10, 12),
  key.title="", key.xlab="", ColSideColors = col, labCol= "")
dev.off()

# Extract cluster information for each case
write.table(cbind(clusters, col), "clusters_Module17.txt", quote = FALSE, row.names = TRUE, col.names = NA, sep = "\t")

# Sum of Z-score
Order <- c2$labels[c2$order]
OK <- cbind(colSums(DATA3)[Order], clusters[Order])
colnames(OK) <- c("Sum_Zscore", "Cluster")
write.table(OK, "Zscore_ClusterOrder_Module17.txt", quote = FALSE, row.names = TRUE, col.names = NA, sep = "\t")

rm(list = ls())
gc();gc()
