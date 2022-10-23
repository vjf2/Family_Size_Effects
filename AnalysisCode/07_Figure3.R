#############
#### Figure 3
#############

# Load standard form model
load("IntermediateData/standardformmod.RData")

# Exponentiate coefficients
exp_coef <- exp(coef(glmod))[2:4]

# Calculate profiled confidence intervals
ci <- confint(glmod, parm = c("maternal", "paternal", "distant"))

# Exponentiate confidence intervals
ci2 <- exp(ci)

windows()
# pdf(file = "TablesFigures/Figure3.pdf", width = 5, height = 4)
par(mar = c(5, 5, 2, 2))

plot(ci2, type = "n",
     ylim = c(0, 9), xlim = c(0.5, 3.5),
     xaxt = "n", yaxt = "n",
     xlab = "Kinship Category",
     ylab = "Odds Ratio (95% CI): Significant Affiliation",
     cex.lab = 1)

abline(h = 1, lty = 5, col = "grey")

segments(1, ci2[1, 1], 1, ci2[1, 2], col = "black")

  segments(0.95, ci2[1, 1], 1.05, ci2[1, 1], col = "black")
  segments(0.95, ci2[1, 2], 1.05, ci2[1, 2], col = "black")

segments(2, ci2[2, 1], 2, ci2[2, 2], col = "black")

  segments(1.95, ci2[2, 1], 2.05, ci2[2, 1], col = "black")
  segments(1.95, ci2[2, 2], 2.05, ci2[2, 2], col = "black")

segments(3, ci2[3, 1], 3, ci2[3, 2], col = "black")

  segments(2.95, ci2[3, 1], 3.05, ci2[3, 1], col = "black")
  segments(2.95, ci2[3, 2], 3.05, ci2[3, 2], col = "black")

points(1:3, exp_coef, pch = 15, cex = 1.5)

axis(1, at = 1:3, labels = c("Matrilineal", "Non-matrilineal", "Distant"))
axis(2, las = 1)

dev.off()
