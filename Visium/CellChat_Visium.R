library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)

library(tidyverse)
library(reshape2)

library(gplots)

set.seed(1234)

packageVersion("Seurat")
packageVersion("CellChat")

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 100 * 1024^3)

# Load a merged Seurat object
lung <- readRDS("lung_visium_merge_harmony_annotation.rds")

lung[["Spatial"]] <- JoinLayers(lung[["Spatial"]])
DefaultAssay(lung) <- "Spatial"
lung

# Please read the sample list (Sample ID, Subtype)
DATA <- read.table("<LIST>", header=T)
case <- DATA$Case

for(i in case) {
      sub <- subset(lung, subset = original == i)
      sub[["samples"]] <- i
      sub$annotation = droplevels(sub$annotation, exclude = setdiff(levels(sub$annotation), unique(sub$annotation)))

      # Set paramters for Visium (https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html)
      dir <- case$DIR[case$ID == i]
      scalefactors = jsonlite::fromJSON(txt = file.path('scalefactors_json.json'))
      spot.size = 65 # use the typical human cell size
      conversion.factor = spot.size / scalefactors$spot_diameter_fullres
      spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

      # Set coordinates
      spatial.locs = GetTissueCoordinates(sub, scale = NULL, cols = c("imagerow", "imagecol"))

      # Creating a CellChat object
      cellchat <- createCellChat(object = sub, assay = "Spatial", group.by = "annotation", datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
      cellchat
      
      # Set the database
      CellChatDB <- CellChatDB.human
      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
      cellchat@DB <- CellChatDB.use
      
      # CellChat
      cellchat <- subsetData(cellchat)
      future::plan("multisession", workers = 4)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
      cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, contact.dependent = TRUE, contact.range = 100)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat) #pathway level
      
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
