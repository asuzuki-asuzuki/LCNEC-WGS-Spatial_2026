library(survival)
library(dplyr)

packageVersion("survival")

# Load survival data with expression patterns of Module 17
data <- read.table("exp_module17_surv.txt", header=T)
data <- data %>%  mutate(Censor = if_else(SurvivalCensor == 0, NA, SurvivalCensor))
data <- data %>%  mutate(Censor = if_else(SurvivalCensor == "alive", 0, 1))

# Survival analysis among 4 clusters
data$DP <- factor(data$DP, levels = c("High", "Middle-High","Middle-Low", "Low"))
KM <- survfit(Surv(SurvivalMonths, Censor) ~ DP, data=data)

table(data$DP[which(!is.na(data$SurvivalMonths))])

# Calculation of p-value (log-rank test)
sd_all <- survdiff(Surv(SurvivalMonths, Censor) ~ DP, data=data)
p_val_all <- 1 - pchisq(sd_all$chisq, length(sd_all$n) - 1)
p_label_all <- format.pval(p_val_all, digits = 3)

# Calculation of p-value each (log-rank test)
DPs <- levels(data$DP)
pairs <- combn(DPs, 2)
results <- list()

for(i in 1:ncol(pairs)){
  g1 <- pairs[1, i]
  g2 <- pairs[2, i]
  sub_data <- subset(data, DP %in% c(g1, g2))
  sd <- survdiff(Surv(SurvivalMonths, Censor) ~ DP, data = sub_data)
  p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)  
  results[[paste(g1, "vs", g2)]] <- p_val
}

# FDR (Benjamini & Hochberg method) correction
p_raw <- unlist(results)
p_adj <- p.adjust(p_raw, method = "fdr")
final_results <- data.frame(Comparison = names(p_raw), Raw_P = p_raw, Adjusted_P = p_adj)
write.table(final_results, file = "p_value_ajusted_module17_exp.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# KM Plot
pdf (file = "module17_exp_KM.pdf", width = 6, height = 6)
plot(KM, lty=c(2, 1, 1, 1), col = c("red", "orange", "lightblue", "blue"), lwd = 2,  xscale=12, xlab="Survival time (years)", ylab="Proportion Surviving", main = "DP-associated module", mark.time=TRUE)
legend ("topright", lty=c(2, 1, 1, 1), legend=levels(data$DP), col = c("red", "orange", "lightblue", "blue"))
legend("bottomleft", legend = c(paste("Log-rank p =", p_label_all)), bty = "n")
dev.off()

