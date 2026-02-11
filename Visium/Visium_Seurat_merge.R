library(Seurat) # v5.0.3
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3)

set.seed(1234)
packageVersion("Seurat")

# Load Seurat objects
t1 <- readRDS("<DIR1>/lung_visium.rds")
t2 <- readRDS("<DIR2>/lung_visium.rds")
t3 <- readRDS("<DIR3>/lung_visium.rds")
t4 <- readRDS("<DIR4>/lung_visium.rds")
t5 <- readRDS("<DIR5>/lung_visium.rds")
t6 <- readRDS("<DIR6>/lung_visium.rds")
t7 <- readRDS("<DIR7>/lung_visium.rds")
t9 <- readRDS("<DIR9>/lung_visium.rds")
t10 <- readRDS("<DIR10>/lung_visium.rds")
t11 <- readRDS("<DIR11>/lung_visium.rds")
t13 <- readRDS("<DIR13>/lung_visium.rds")
t14 <- readRDS("<DIR14>/lung_visium.rds")
j10 <- readRDS("<DIR15>/lung_visium.rds")
j23 <- readRDS("<DIR16>/lung_visium.rds")
j25 <- readRDS("<DIR17>/lung_visium.rds")
j27 <- readRDS("<DIR18>/lung_visium.rds")
j28 <- readRDS("<DIR19>/lung_visium.rds")
j29 <- readRDS("<DIR20>/lung_visium.rds")
j33 <- readRDS("<DIR21>/lung_visium.rds")
# Total of 19 objects

t1$original <- "<SAMPLE1>"
t2$original <- "<SAMPLE2>"
t3$original <- "<SAMPLE3>"
t4$original <- "<SAMPLE4>"
t5$original <- "<SAMPLE5>"
t6$original <- "<SAMPLE6>"
t7$original <- "<SAMPLE7>"
t9$original <- "<SAMPLE9>"
t10$original <- "<SAMPLE10>"
t11$original <- "<SAMPLE11>"
t13$original <- "<SAMPLE13>"
t14$original <- "<SAMPLE14>"
j10$original <- "<SAMPLE15>"
j23$original <- "<SAMPLE16>"
j25$original <- "<SAMPLE17>"
j27$original <- "<SAMPLE18>"
j28$original <- "<SAMPLE19>"
j29$original <- "<SAMPLE20>"
j33$original <- "<SAMPLE21>"
# Total of 19 sections

# Merge
lung.merge <- merge(t1, y = c(t2, t3, t4, t5, t6, t7, t9, t10, t11, t13, t14, j10, j23, j25, j27, j28, j29, j33), add.cell.ids = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t9", "t10", "t11", "t13", "t14", "j10", "j23", "j25", "j27", "j28", "j29", "j33"), project = "LCNEC")

# Number of spots
length(lung.merge@meta.data$original)
table(lung.merge@meta.data$original)

# Save the object
saveRDS(lung.merge, file = "lung_visium_merge.rds")

rm(list = ls())
gc();gc()
