#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-04
# Description: Relative importance decomposition of ILS, GTEE, and IH
#              on gene tree variation (gCF).
#              Includes LOESS smoothing with 95% CI for ILS and GTEE.

# ====== Main Script ======

library(relaimpo)
library(ggplot2)
library(gridExtra)

args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript relaimpo.R <input.csv> [suffix]")
}

# 第二个参数为可选后缀，用于输出文件名
if (length(args) >= 2) {
  suffix <- paste0("_", args[2])
} else {
  suffix <- ""
}

x <- read.csv(args[1], header = TRUE, sep = ",")

# 移除含 NA 的行
data <- x[, c('genetree_var', 'ILS', 'GTEE', 'IH')]
data <- data[complete.cases(data), ]
cat(sprintf("有效节点数: %d\n", nrow(data)))

# ====== 1. Relative Importance (relaimpo) ======

cat("\n===== Relative Importance Decomposition =====\n")
relimp <- calc.relimp(data, type = "lmg", rela = TRUE)
print(relimp)

# Bootstrap for confidence intervals
set.seed(42)
bootimpo.result <- boot.relimp(
  genetree_var ~ ILS + GTEE + IH,
  data = data, b = 1000,
  type = c("lmg", "last", "first", "pratt"),
  rank = TRUE, diff = TRUE, rela = TRUE
)

cat("\n===== Bootstrap Results =====\n")
print(booteval.relimp(bootimpo.result, lev = 0.9, nodiff = TRUE))

# Plot relaimpo
pdf(paste0("relimpo_bootstrap", suffix, ".pdf"), width = 10, height = 8)
plot(booteval.relimp(bootimpo.result))
dev.off()
cat(sprintf("已保存: relimpo_bootstrap%s.pdf\n", suffix))

# ====== 2. Linear Model ======

cat("\n===== Linear Model =====\n")
linmod_basic <- lm(genetree_var ~ ILS + GTEE + IH, data = data)
print(summary(linmod_basic))

# ====== 3. Pearson Correlation ======

cat("\n===== Pearson Correlations =====\n")
cor_ils <- cor.test(data$ILS, data$genetree_var, method = "pearson")
cor_gtee <- cor.test(data$GTEE, data$genetree_var, method = "pearson")
cor_ih <- cor.test(data$IH, data$genetree_var, method = "pearson")

cat(sprintf("  ILS  vs genetree_var: r = %.4f, p = %.2e\n",
            cor_ils$estimate, cor_ils$p.value))
cat(sprintf("  GTEE vs genetree_var: r = %.4f, p = %.2e\n",
            cor_gtee$estimate, cor_gtee$p.value))
cat(sprintf("  IH   vs genetree_var: r = %.4f, p = %.2e\n",
            cor_ih$estimate, cor_ih$p.value))

# ====== 4. LOESS Plots with 95% CI ======

# ILS vs genetree_var
p_ils <- ggplot(data, aes(x = ILS, y = genetree_var)) +
  geom_point(color = "black", size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", fill = "grey70",
              alpha = 0.3, se = TRUE, level = 0.95) +
#   annotate("text",
#            x = min(data$ILS) + 0.01 * diff(range(data$ILS)),
#            y = max(data$genetree_var),
#            label = sprintf("r = %.3f, p = %.2e", cor_ils$estimate, cor_ils$p.value),
#            hjust = 0, size = 4, color = "darkred") +
    annotate("text",
         x = -Inf, y = Inf,
         label = sprintf("r = %.3f, p = %.2e", cor_ils$estimate, cor_ils$p.value),
         hjust = -0.1, vjust = 1.5,
         size = 4, color = "darkred") +
  labs(x = "ILS (theta)",
       y = "Gene tree variation (100 - gCF)",
       title = "ILS vs Gene Tree Variation") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# GTEE vs genetree_var
p_gtee <- ggplot(data, aes(x = GTEE, y = genetree_var)) +
  geom_point(color = "black", size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", fill = "grey70",
              alpha = 0.3, se = TRUE, level = 0.95) +
#   annotate("text",
#            x = min(data$GTEE) + 0.01 * diff(range(data$GTEE)),
#            y = max(data$genetree_var),
#            label = sprintf("r = %.3f, p = %.2e", cor_gtee$estimate, cor_gtee$p.value),
#            hjust = 0, size = 4, color = "darkred") +
    annotate("text",
         x = -Inf, y = Inf,
         label = sprintf("r = %.3f, p = %.2e", cor_gtee$estimate, cor_gtee$p.value),
         hjust = -0.1, vjust = 1.5,
         size = 4, color = "darkred") +
  labs(x = "GTEE (100 - recovery %)",
       y = "Gene tree variation (100 - gCF)",
       title = "GTEE vs Gene Tree Variation") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# IH vs genetree_var
p_ih <- ggplot(data, aes(x = IH, y = genetree_var)) +
  geom_point(color = "black", size = 1.5, alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", fill = "grey70",
              alpha = 0.3, se = TRUE, level = 0.95) +
#   annotate("text",
#            x = min(data$IH) + 0.01 * diff(range(data$IH)),
#            y = max(data$genetree_var),
#            label = sprintf("r = %.3f, p = %.2e", cor_ih$estimate, cor_ih$p.value),
#            hjust = 0, size = 4, color = "darkred") +
  annotate("text",
         x = -Inf, y = Inf,
         label = sprintf("r = %.3f, p = %.2e", cor_ih$estimate, cor_ih$p.value),
         hjust = -0.1, vjust = 1.5,
         size = 4, color = "darkred") +
  labs(x = "IH (reticulation index)",
       y = "Gene tree variation (100 - gCF)",
       title = "IH vs Gene Tree Variation") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# 保存 LOESS 图
pdf(paste0("relimpo_loess_correlation", suffix, ".pdf"), width = 18, height = 5)
grid.arrange(p_ils, p_gtee, p_ih, ncol = 3)
dev.off()
cat(sprintf("已保存: relimpo_loess_correlation%s.pdf\n", suffix))

# ====== 5. Summary ======

cat("\n===== Summary =====\n")
cat(sprintf("R-squared (linear model): %.4f\n", summary(linmod_basic)$r.squared))
cat("\nRelative contributions (LMG, %):\n")
cat(sprintf("  ILS:  %.1f%%\n", relimp@lmg[1] * 100))
cat(sprintf("  GTEE: %.1f%%\n", relimp@lmg[2] * 100))
cat(sprintf("  IH:   %.1f%%\n", relimp@lmg[3] * 100))
cat("\nDone.\n")
