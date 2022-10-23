#############
#### Figure 4
#############

library(viridisLite)

load("IntermediateData/model_results.RData") # generated from 05_Gregariousness.R

# create dataframe for model predictions
newdata <- with(msm, expand.grid(avail_assoc = mean(avail_assoc),
                                 matkins = 0:4,
                                 patkins = mean(patkins),
                                 distkins = mean(distkins),
                                 ForagingType = as.factor(c(1:4)),
                                 provisioned = "N",
                                 adult_age = mean(adult_age),
                                 sightings = mean(sightings),
                                 adult_yrs = mean(adult_yrs)
))

ilink <- family(tmod)$linkinv
p <- predict(tmod, newdata, se.fit = TRUE)

t_predict <- cbind(newdata,
                   res  = ilink(p$fit),
                   Ures = ilink(p$fit + (1.96 * p$se.fit)),
                   Lres = ilink(p$fit - (1.96 * p$se.fit)))

p2 <- predict(emod, newdata, se.fit = TRUE)

e_predict <- cbind(newdata,
                   res  = p2$fit,
                   Ures = p2$fit + (1.96 * p2$se.fit),
                   Lres = p2$fit - (1.96 * p2$se.fit))

p3 <- predict(umod, newdata, se.fit = TRUE)

u_predict <- cbind(newdata,
                   res  = ilink(p3$fit),
                   Ures = ilink(p3$fit + (1.96 * p3$se.fit)),
                   Lres = ilink(p3$fit - (1.96 * p3$se.fit)))

p4 <- predict(cmod, newdata, se.fit = TRUE)

c_predict <- cbind(newdata,
                   res  = p4$fit,
                   Ures = p4$fit + (1.96 * p4$se.fit),
                   Lres = p4$fit - (1.96 * p4$se.fit))

set.seed(286567440) # to preserve jitter

msm$jmatkins <- jitter(msm$matkins, factor = 0.75)
msm$jmatkins <- ifelse(msm$jmatkins <= 0, 0, msm$jmatkins)
msm$jmatkins <- ifelse(msm$jmatkins >= 4, 4, msm$jmatkins)

msm$forcolor[msm$ForagingType == 1] <- viridis(4)[1] # purple viridis(4)[1]
msm$forcolor[msm$ForagingType == 2] <- viridis(4)[2] # blue viridis(4)[2]
msm$forcolor[msm$ForagingType == 3] <- viridis(4)[4] # yellow viridis(4)[4]
msm$forcolor[msm$ForagingType == 4] <- viridis(4)[3] # green viridis(4)[3]

msm$shape <- ifelse(msm$provisioned == "Y", 17, 16)

# plot

windows()
# pdf("TablesFigures/Figure4.pdf", width = 6)
par(mar = c(2, 4.1, 0.8, 1), oma = c(4, 2, 0, 2), xpd = NA)

mat <- matrix(rep(c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10), each = 4), ncol = 2)

mat <- mat[-c(10:12), ]

layout(mat)


plot(msm$jmatkins, msm$tot_affil, ylim = c(0, 18),
     xlab = NA,
     ylab = "Total Affiliations",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor("black", alpha.f = 0.7),
     cex = 1.25)

axis(1)
axis(2, las = 1)

# making the predictions
p <- aggregate(res ~ matkins, data = t_predict, mean)

l <- aggregate(Lres ~ matkins, data = t_predict, mean)
lines(l$matkins, l$Lres, lty = 2)

u <- aggregate(Ures ~ matkins, data = t_predict, mean)
lines(u$matkins, u$Ures, lty = 2)

polygon(c(l$matkins, rev(l$matkins)), c(l$Lres, rev(u$Ures)),
        col = adjustcolor("grey", alpha.f = 0.6), border = NA)

lines(p$matkins, p$res, lty = 1)

# unrelated affiliations

plot(msm$jmatkins, msm$unrel_affil, ylim = c(0, 18),
     xlab = NA,
     ylab = "Unrelated Affiliations",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor("black", alpha.f = 0.7),
     cex = 1.25)

legend("topright", pch = c(1, 2), legend = c("Non-provisioned", "Provisioned"), bty = "n", cex = 0.9)

axis(1)
axis(2, las = 1)
text(0.5, 17, "N.S.", font = 2)

# eigenvector centrality
plot(msm$jmatkins, scale(msm$resid_eigen), ylim = c(-3.3, 3.3),
     xlab = NA,
     ylab = "Eigenvector Centrality",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor("black", alpha.f = 0.7),
     cex = 1.25)

axis(1)
axis(2, las = 1)

# model predictions
p2 <- aggregate(res ~ matkins, data = e_predict, mean)

l <- aggregate(Lres ~ matkins, data = e_predict, mean)
lines(l$matkins, l$Lres, lty = 2)

u <- aggregate(Ures ~ matkins, data = e_predict, mean)
lines(u$matkins, u$Ures, lty = 2)

polygon(c(l$matkins, rev(l$matkins)), c(l$Lres, rev(u$Ures)),
        col = adjustcolor("grey", alpha.f = 0.6), border = NA)

lines(p2$matkins, p2$res, lty = 1)

# closeness centrality
plot(msm$jmatkins, scale(msm$resid_closeness), ylim = c(-3.3, 3.3),
     xlab = NA,
     ylab = "Closeness Centrality",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor("black", alpha.f = 0.7),
     cex = 1.25)

axis(1)
axis(2, las = 1)
text(0.5, 3, "N.S.", font = 2)

## Two blank plots for labels
plot.new()
plot.new()

## Foraging Type plots start here

# create dataframe for model predictions
newdata <- with(msm, expand.grid(avail_assoc = mean(avail_assoc),
                                 matkins = mean(matkins),
                                 patkins = mean(patkins),
                                 distkins = mean(distkins),
                                 ForagingType = as.factor(c(1:4)),
                                 provisioned = "N",
                                 adult_age = mean(adult_age),
                                 sightings = mean(sightings),
                                 adult_yrs = mean(adult_yrs)
))

ilink <- family(tmod)$linkinv
p <- predict(tmod, newdata, se.fit = TRUE)

t_predict <- cbind(newdata,
                   res  = ilink(p$fit),
                   Ures = ilink(p$fit + (1.96 * p$se.fit)),
                   Lres = ilink(p$fit - (1.96 * p$se.fit)))

p2 <- predict(emod, newdata, se.fit = TRUE)

e_predict <- cbind(newdata,
                   res  = p2$fit,
                   Ures = p2$fit + (1.96 * p2$se.fit),
                   Lres = p2$fit - (1.96 * p2$se.fit))

p3 <- predict(umod, newdata, se.fit = TRUE)

u_predict <- cbind(newdata,
                   res  = ilink(p3$fit),
                   Ures = ilink(p3$fit + (1.96 * p3$se.fit)),
                   Lres = ilink(p3$fit - (1.96 * p3$se.fit)))

p4 <- predict(cmod, newdata, se.fit = TRUE)

c_predict <- cbind(newdata,
                   res  = p4$fit,
                   Ures = p4$fit + (1.96 * p4$se.fit),
                   Lres = p4$fit - (1.96 * p4$se.fit))

msm$jForagingType <- jitter(as.numeric(msm$ForagingType), factor = 0.75)

plot(msm$jForagingType, msm$tot_affil, ylim = c(0, 18),
     xlab = NA,
     ylab = "Total Affiliations",
     xaxt = "n",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor(msm$forcolor, alpha.f = 0.5),
     cex = 1.25)

points(msm$jForagingType[which(msm$provisioned == "Y")],
       msm$tot_affil[which(msm$provisioned == "Y")],
       pch = 2, col = adjustcolor("black", alpha.f = 0.2))

axis(1, at = 1:4, labels = c("TDPD", "MIXED", "SPONGE", "SEAGRASS"), cex.axis = 0.9)
axis(2, las = 1)

# making the predictions
p <- aggregate(res ~ ForagingType, data = t_predict, mean)

l <- aggregate(Lres ~ ForagingType, data = t_predict, mean)

u <- aggregate(Ures ~ ForagingType, data = t_predict, mean)

points(1:4, p$res, pch = 16, cex = 1.5)
segments(1:4, l$Lres, 1:4, u$Ures, lwd = 1.5)

# unrelated affiliations

plot(msm$jForagingType, msm$unrel_affil, ylim = c(0, 18),
     xlab = NA,
     ylab = "Unrelated Affiliations",
     xaxt = "n",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor(msm$forcolor, alpha.f = 0.5),
     cex = 1.25)

points(msm$jForagingType[which(msm$provisioned == "Y")],
       msm$unrel_affil[which(msm$provisioned == "Y")],
       pch = 2, col = adjustcolor("black", alpha.f = 0.2))

legend("topright", pch = c(1, 2), legend = c("Non-provisioned", "Provisioned"), bty = "n", cex = 0.9)

# model predictions
p2 <- aggregate(res ~ ForagingType, data = u_predict, mean)

l <- aggregate(Lres ~ ForagingType, data = u_predict, mean)

u <- aggregate(Ures ~ ForagingType, data = u_predict, mean)

axis(1, at = 1:4, labels = c("TDPD", "MIXED", "SPONGE", "SEAGRASS"), cex.axis = 0.9)
points(1:4, p2$res, pch = 16, cex = 1.5)
segments(1:4, l$Lres, 1:4, u$Ures, lwd = 1.5)

axis(2, las = 1)

# eigenvector centrality
plot(msm$jForagingType, scale(msm$resid_eigen), ylim = c(-3.3, 3.3),
     xlab = NA,
     ylab = "Eigenvector Centrality",
     xaxt = "n",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor(msm$forcolor, alpha.f = 0.5),
     cex = 1.25)

points(msm$jForagingType[which(msm$provisioned == "Y")],
       scale(msm$resid_eigen)[which(msm$provisioned == "Y")],
       pch = 2, col = adjustcolor("black", alpha.f = 0.2))

axis(2, las = 1)

# model predictions
p2 <- aggregate(res ~ ForagingType, data = e_predict, mean)

l <- aggregate(Lres ~ ForagingType, data = e_predict, mean)

u <- aggregate(Ures ~ ForagingType, data = e_predict, mean)

axis(1, at = 1:4, labels = c("TDPD", "MIXED", "SPONGE", "SEAGRASS"), cex.axis = 0.9)
points(1:4, p2$res, pch = 16, cex = 1.5)
segments(1:4, l$Lres, 1:4, u$Ures, lwd = 1.5)

# closeness centrality
plot(msm$jForagingType, scale(msm$resid_closeness), ylim = c(-3.3, 3.3),
     xlab = NA,
     ylab = "Closeness Centrality",
     xaxt = "n",
     yaxt = "n",
     pch = msm$shape,
     col = adjustcolor(msm$forcolor, alpha.f = 0.5),
     cex = 1.25)

points(msm$jForagingType[which(msm$provisioned == "Y")],
       scale(msm$resid_closeness)[which(msm$provisioned == "Y")],
       pch = 2, col = adjustcolor("black", alpha.f = 0.2))

axis(2, las = 1)


# model predictions
p2 <- aggregate(res ~ ForagingType, data = c_predict, mean)

l <- aggregate(Lres ~ ForagingType, data = c_predict, mean)

u <- aggregate(Ures ~ ForagingType, data = c_predict, mean)

axis(1, at = 1:4, labels = c("TDPD", "MIXED", "SPONGE", "SEAGRASS"), cex.axis = 0.9)

points(1:4, p2$res, pch = 16, cex = 1.5)
segments(1:4, l$Lres, 1:4, u$Ures, lwd = 1.5)

text(0, 15, "Number of Matrilineal Close Kin", cex = 1.25)
text(0, -6.5, "Foraging Type", cex = 1.25)

text(-5, 33.9, "(a)", cex = 1.5)
text(-5, 12.75, "(b)", cex = 1.5)

dev.off()
