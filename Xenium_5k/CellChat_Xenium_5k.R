library(Seurat)  # v5.4.0
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gplots)
library(CellChat)
library(Matrix)
library(patchwork)
options(stringsAsFactors = F)

set.seed(1234)

packageVersion("Seurat")
packageVersion("CellChat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a merged Seurat object
lung <- readRDS("lung_xenium_merge_sketch_ProjectData_annotation.rds")

lung[["Xenium"]] <- JoinLayers(lung[["Xenium"]])
DefaultAssay(lung) <- "Xenium"
lung

# Please read the sample list (Sample ID, Subtype)
DATA <- read.table("<LIST>", header=T)
case <- DATA$Case

for(i in case) {
	  # Subset the object in each section
      sub <- subset(lung, subset = original == i)
      sub[["samples"]] <- i
      sub$merge_annotation = droplevels(sub$merge_annotation, exclude = setdiff(levels(sub$merge_annotation), unique(sub$merge_annotation)))

      # Set parameters for Xenium (https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html)
      conversion.factor = 1
      spot.size = 10
      spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

      # Set coordinates
      coords <- GetTissueCoordinates(sub)
      coords <- coords[, c("x", "y")]
      rownames(coords) <- colnames(sub)

      # Creating a CellChat object
      data.input = LayerData(sub, assay = "Xenium", layer = "data")
      meta <- as.data.frame(sub@meta.data)
      cellchat <- createCellChat(object = as.sparse(data.input), meta = meta, group.by = "merge_annotation", datatype = "spatial", coordinates = coords, spatial.factors = spatial.factors)
      cellchat
      
      # Set the database
      CellChatDB <- CellChatDB.human
      CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), non_protein = F)
      cellchat@DB <- CellChatDB.use
      
      # CellChat
      cellchat <- subsetData(cellchat)
      future::plan("multisession", workers = 7)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
      cellchat <- computeCommunProb(cellchat, distance.use = TRUE, type = "truncatedMean", trim = 0.1, interaction.range = 250, scale.distance = 1, contact.dependent = TRUE, contact.range = 10)
      cellchat <- filterCommunication(cellchat, min.cells = 500)
      cellchat <- computeCommunProbPathway(cellchat) # pathway level
      
      # Output to a data frame
      df <- subsetCommunication(cellchat)
      write.csv(df, paste0("cellchat/", i, "_cellchat_each.csv"), row.names = FALSE)
      df <- subsetCommunication(cellchat, slot.name = "netP")
      write.csv(df, paste0("cellchat/", i, "_cellchat_pathway_each.csv"), row.names = FALSE)

      # Save a CellChat object
      saveRDS(cellchat, file=paste0("cellchat/", i, "_cellchat.rds"))
}

rm(list = ls())
gc();gc()
