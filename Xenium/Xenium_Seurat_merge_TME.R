library(Seurat) # v5.0.3
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3) # 100 Gb

set.seed(1234)
packageVersion("Seurat")

# Load Seurat objects
lung_t1 <- readRDS("<DIR1>/lung_niche_TME.rds")
lung_t2 <- readRDS("<DIR2>/lung_niche_TME.rds")
lung_t3 <- readRDS("<DIR3>/lung_niche_TME.rds")
lung_t4 <- readRDS("<DIR4>/lung_niche_TME.rds")
lung_t5 <- readRDS("<DIR5>/lung_niche_TME.rds")
lung_t6 <- readRDS("<DIR6>/lung_niche_TME.rds")
lung_t7 <- readRDS("<DIR7>/lung_niche_TME.rds")
lung_t8 <- readRDS("<DIR8>/lung_niche_TME.rds")
lung_t9 <- readRDS("<DIR9>/lung_niche_TME.rds")
lung_t10 <- readRDS("<DIR10>/lung_niche_TME.rds")
lung_t11 <- readRDS("<DIR11>/lung_niche_TME.rds")
lung_t12 <- readRDS("<DIR12>/lung_niche_TME.rds")
lung_t13 <- readRDS("<DIR13>/lung_niche_TME.rds")
lung_t14 <- readRDS("<DIR14>/lung_niche_TME.rds")
lung_t15 <- readRDS("<DIR15>/lung_niche_TME.rds")
lung_t16_a <- readRDS("<DIR16a>/lung_niche_TME.rds")
lung_t16_b <- readRDS("<DIR16b>/lung_niche_TME.rds")
lung_t17 <- readRDS("<DIR17>/lung_niche_TME.rds")
lung_t18 <- readRDS("<DIR18>/lung_niche_TME.rds")
lung_t19_a <- readRDS("<DIR19a>/lung_niche_TME.rds")
lung_t19_b <- readRDS("<DIR19b>/lung_niche_TME.rds")
lung_t20 <- readRDS("<DIR20>/lung_niche_TME.rds")
lung_t21 <- readRDS("<DIR21>/lung_niche_TME.rds")
lung_t22 <- readRDS("<DIR22>/lung_niche_TME.rds")
lung_t23 <- readRDS("<DIR23>/lung_niche_TME.rds")
lung_t24 <- readRDS("<DIR24>/lung_niche_TME.rds")
lung_t25 <- readRDS("<DIR25>/lung_niche_TME.rds")
lung_t25_2 <- readRDS("<DIR25_2>/lung_niche_TME.rds")
lung_t26 <- readRDS("<DIR26>/lung_niche_TME.rds")
lung_t27 <- readRDS("<DIR27>/lung_niche_TME.rds")
lung_j10 <- readRDS("<DIR28>/lung_niche_TME.rds")
lung_j23 <- readRDS("<DIR29>/lung_niche_TME.rds")
lung_j25 <- readRDS("<DIR30>/lung_niche_TME.rds")
lung_j27 <- readRDS("<DIR31>/lung_niche_TME.rds")
lung_j28 <- readRDS("<DIR32>/lung_niche_TME.rds")
lung_j29 <- readRDS("<DIR33>/lung_niche_TME.rds")
lung_j33 <- readRDS("<DIR34>/lung_niche_TME.rds")
# Total of 37 objects

# Input sample names to "original" in each object
lung_t1$original <- "<SAMPLE1>"
lung_t2$original <- "<SAMPLE2>"
lung_t3$original <- "<SAMPLE3>"
lung_t4$original <- "<SAMPLE4>"
lung_t5$original <- "<SAMPLE5>"
lung_t6$original <- "<SAMPLE6>"
lung_t7$original <- "<SAMPLE7>"
lung_t8$original <- "<SAMPLE8>"
lung_t9$original <- "<SAMPLE9>"
lung_t10$original <- "<SAMPLE10>"
lung_t11$original <- "<SAMPLE11>"
lung_t12$original <- "<SAMPLE12>"
lung_t13$original <- "<SAMPLE13>"
lung_t14$original <- "<SAMPLE14>"
lung_t15$original <- "<SAMPLE15>"
lung_t16_a$original <- "<SAMPLE16a>"
lung_t16_b$original <- "<SAMPLE16b>"
lung_t17$original <- "<SAMPLE17>"
lung_t18$original <- "<SAMPLE18>"
lung_t19_a$original <- "<SAMPLE19a>"
lung_t19_b$original <- "<SAMPLE19b>"
lung_t20$original <- "<SAMPLE20>"
lung_t21$original <- "<SAMPLE21>"
lung_t22$original <- "<SAMPLE22>"
lung_t23$original <- "<SAMPLE23>"
lung_t24$original <- "<SAMPLE24>"
lung_t25$original <- "<SAMPLE25>"
lung_t25_2$original <- "<SAMPLE25_2>"
lung_t26$original <- "<SAMPLE26>"
lung_t27$original <- "<SAMPLE27>"
lung_j10$original <- "<SAMPLE28>"
lung_j23$original <- "<SAMPLE29>"
lung_j25$original <- "<SAMPLE30>"
lung_j27$original <- "<SAMPLE31>"
lung_j28$original <- "<SAMPLE32>"
lung_j29$original <- "<SAMPLE33>"
lung_j33$original <- "<SAMPLE34>"
# Total of 37 sections

# Merge
lung.merge <- merge(lung_t1, y = c(lung_t2, lung_t3, lung_t4, lung_t5, lung_t6, lung_t7, lung_t8, lung_t9, lung_t10, lung_t11, lung_t12, lung_t13, lung_t14, lung_t15, lung_t16_a, lung_t16_b, lung_t17, lung_t18, lung_t19_a, lung_t19_b, lung_t20, lung_t21, lung_t22, lung_t23, lung_t24, lung_t25, lung_t25_2, lung_t26, lung_t27, lung_j10, lung_j23, lung_j25, lung_j27, lung_j28, lung_j29, lung_j33), add.cell.ids = c("lung_t1", "lung_t2", "lung_t3", "lung_t4", "lung_t5", "lung_t6", "lung_t7", "lung_t8", "lung_t9", "lung_t10", "lung_t11", "lung_t12", "lung_t13", "lung_t14", "lung_t15", "lung_t16_a", "lung_t16_b", "lung_t17", "lung_t18", "lung_t19_a", "lung_t19_b", "lung_t20", "lung_t21", "lung_t22", "lung_t23", "lung_t24", "lung_t25", "lung_t25_2", "lung_t26", "lung_t27", "lung_j10", "lung_j23", "lung_j25", "lung_j27", "lung_j28", "lung_j29", "lung_j33"), project = "LCNEC")

# Number of cells
length(lung.merge@meta.data$original)
table(lung.merge@meta.data$original)

# Save the object
saveRDS(lung.merge, file = "lung_niche_TME_merge.rds")

rm(list = ls())
gc();gc()
