library(survival)
library(dplyr)

# Load the survival data with mutation status of ADAMTS family genes
data <- read.table("mut_ADAMTS.txt", header=T)

data <- data %>%  mutate(Censor = if_else(SurvivalCensor == 0, NA, SurvivalCensor))
data <- data %>%  mutate(Censor = if_else(SurvivalCensor == "alive", 0, 1))

# Survival analysis in ADAMTS20
data$ADAMTS20 <- factor(data$ADAMTS20, levels = c(1, 0), labels=c("mutation(+)", "mutation(-)"))
KM <- survfit(Surv(SurvivalMonths, Censor) ~ ADAMTS20, data=data)
sd <- survdiff(Surv(SurvivalMonths, Censor) ~ ADAMTS20, data=data)
p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_label <- format.pval(p_val, digits = 3)

pdf (file = "mut_KM_ADAMTS20.pdf", width = 5, height = 5)
plot(KM, lty=2:1, col = c("red", "black"), lwd = 2,  xscale=12, xlab="Survival time (years)", ylab="Proportion Surviving", main = "ADAMTS20", mark.time=TRUE)
legend ("topright", lty=2:1, legend=levels(data$ADAMTS20), col = c("red", "black"))
legend("bottomleft", legend = c(paste("Log-rank p =", p_label)), bty = "n")
dev.off()

# Survival analysis in ADAMTS2
data$ADAMTS2 <- factor(data$ADAMTS2, levels = c(1, 0), labels=c("mutation(+)", "mutation(-)"))
KM <- survfit(Surv(SurvivalMonths, Censor) ~ ADAMTS2, data=data)
sd <- survdiff(Surv(SurvivalMonths, Censor) ~ ADAMTS2, data=data)
p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_label <- format.pval(p_val, digits = 3)

pdf (file = "mut_KM_ADAMTS2.pdf", width = 5, height = 5)
plot(KM, lty=2:1, col = c("red", "black"), lwd = 2,  xscale=12, xlab="Survival time (years)", ylab="Proportion Surviving", main = "ADAMTS2", mark.time=TRUE)
legend ("topright", lty=2:1, legend=levels(data$ADAMTS2), col = c("red", "black"))
legend("bottomleft", legend = c(paste("Log-rank p =", p_label)), bty = "n")
dev.off()

# Survival analysis in ADAMTS12
data$ADAMTS12 <- factor(data$ADAMTS12, levels = c(1, 0), labels=c("mutation(+)", "mutation(-)"))
KM <- survfit(Surv(SurvivalMonths, Censor) ~ ADAMTS12, data=data)
sd <- survdiff(Surv(SurvivalMonths, Censor) ~ ADAMTS12, data=data)
p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_label <- format.pval(p_val, digits = 3)

pdf (file = "mut_KM_ADAMTS12.pdf", width = 5, height = 5)
plot(KM, lty=2:1, col = c("red", "black"), lwd = 2,  xscale=12, xlab="Survival time (years)", ylab="Proportion Surviving", main = "ADAMTS12", mark.time=TRUE)
legend ("topright", lty=2:1, legend=levels(data$ADAMTS12), col = c("red", "black"))
legend("bottomleft", legend = c(paste("Log-rank p =", p_label)), bty = "n")
dev.off()


# Survival analysis in ADAMTS9
data$ADAMTS9 <- factor(data$ADAMTS9, levels = c(1, 0), labels=c("mutation(+)", "mutation(-)"))
KM <- survfit(Surv(SurvivalMonths, Censor) ~ ADAMTS9, data=data)
sd <- survdiff(Surv(SurvivalMonths, Censor) ~ ADAMTS9, data=data)
p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_label <- format.pval(p_val, digits = 3)

pdf (file = "mut_KM_ADAMTS9.pdf", width = 5, height = 5)
plot(KM, lty=2:1, col = c("red", "black"), lwd = 2,  xscale=12, xlab="Survival time (years)", ylab="Proportion Surviving", main = "ADAMTS9", mark.time=TRUE)
legend ("topright", lty=2:1, legend=levels(data$ADAMTS9), col = c("red", "black"))
legend("bottomleft", legend = c(paste("Log-rank p =", p_label)), bty = "n")
dev.off()
