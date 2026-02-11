library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a Seurat object
lung.merge <- readRDS("lung_sub_sketch_ProjectData_subannotation.rds")

# Subset tumor cell sub-clusters
lung.merge <- subset(lung.merge, subset = subannotation %in% c("NE-1", "NE-2", "NE-3", "NE-4", "NE-5", "NE-6", "NE-7", "NE-8", "NE-9", "Non-NE-p-1", "Non-NE-p-2", "Non-NE-p-3", "Non-NE-p-4", "Non-NE-p-5", "Non-NE-y-1", "Non-NE-y-2", "Non-NE-y-3", "NE/Non-NE-p-1", "NE/Non-NE-y-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3"))

# Xenium assay
DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# Join layers
lung.merge[["Xenium"]] <- JoinLayers(lung.merge[["Xenium"]])

# Set thresholds of molecules for positive
num <- 1
num2 <- 1
num3 <- 1

lung.merge[["TFsig"]] <- NA
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] >= num & lung.merge[["Xenium"]]$counts["ASCL1",] >= num2 & lung.merge[["Xenium"]]$counts["ASCL2",] >= num3)] <- "TP"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] >= num & lung.merge[["Xenium"]]$counts["ASCL1",] >= num2 & lung.merge[["Xenium"]]$counts["ASCL2",] < num3)] <- "DP_YAP1_ASCL1"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] >= num & lung.merge[["Xenium"]]$counts["ASCL1",] < num2 & lung.merge[["Xenium"]]$counts["ASCL2",] >= num3)] <- "DP_YAP1_ASCL2"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] < num & lung.merge[["Xenium"]]$counts["ASCL1",] >= num2 & lung.merge[["Xenium"]]$counts["ASCL2",] >= num3)] <- "DP_ASCL1_ASCL2"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] >= num & lung.merge[["Xenium"]]$counts["ASCL1",] < num2 & lung.merge[["Xenium"]]$counts["ASCL2",] < num3)] <- "YAP1"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] < num & lung.merge[["Xenium"]]$counts["ASCL1",] >= num2 & lung.merge[["Xenium"]]$counts["ASCL2",] < num3)] <- "ASCL1"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] < num & lung.merge[["Xenium"]]$counts["ASCL1",] < num2 & lung.merge[["Xenium"]]$counts["ASCL2",] >= num3)] <- "ASCL2"
lung.merge$TFsig[which(lung.merge[["Xenium"]]$counts["YAP1",] < num & lung.merge[["Xenium"]]$counts["ASCL1",] < num2 & lung.merge[["Xenium"]]$counts["ASCL2",] < num3)] <- "TN"

# Output the count table
OK <- table(lung.merge$original, lung.merge$TFsig)
write.table(OK, file = "count_ASCL1_ASCL2_YAP1.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


##########
#PCA plot#
##########
DATA <- read.table("count_ASCL1_ASCL2_YAP1.txt", header=T)
DATA2 <- DATA/rowSums(DATA)

# Set colors in subtypes
OK <- read.table("id_subtype_NE_nonNE_37.txt", check.names=FALSE, header=T) # please prepare the list of sample names and subtypes (ID: sample names; NE_nonNE: subtypes including "NE", "nonNE", "Mixed" and "DP")
OK$COL <- "NA"
OK$COL[which(OK$NE_nonNE == "NE")] <- "blue" 
OK$COL[which(OK$NE_nonNE == "nonNE")] <- "red"
OK$COL[which(OK$NE_nonNE == "Mixed")] <- "midnightblue"
OK$COL[which(OK$NE_nonNE == "DP")] <- "purple"

# PCA
pca <- prcomp (DATA2, scale.=T)
pcScores <- pca$x[,1:2]

pdf("pca_NE_nonNE.pdf", width=5, height=5)
p <- ggplot(pcScores, aes(PC1, PC2, label = OK$ID)) + geom_point(color = OK$COL, size = 3) + geom_text_repel(size=2.5, max.overlaps=Inf)
p
dev.off()

rm(list = ls())
gc();gc()
