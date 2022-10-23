####################################
#### Model Individual Gregariousness
####################################

library(MASS)

load("IntermediateData/centrality_networks.RData")

# read in focaldata generated from 03_Affiliations.R
affildata <- read.csv("IntermediateData/affiliations.csv", row.names = 1)

min(affildata$HRO[which(affildata$affiliation975 == TRUE)])

# count up total and unrelated affils per female and total possible affiliates
focal_adults$avail_assoc <- NA
focal_adults$tot_affil <- NA
focal_adults$unrel_affil <- NA

for (i in 1:nrow(focal_adults)) {

  ego <- focal_adults[i, "Dolphin.ID"]
  theserows <- affildata[which(affildata$ID1 == ego & affildata$adult_yrsoverlap >= 1), ]
  focal_adults[i, "avail_assoc"] <- sum(theserows$HRO >= 0.25)
  focal_adults[i, "tot_affil"] <- sum(theserows$affiliation975)
  focal_adults[i, "unrel_affil"] <- sum(theserows$affiliation975[theserows$kintype == "unrelated"])
}

# model social data
msm <- focal_adults
msm$sightings <- log(msm$num_sightings_adult)

# set most common foraging type (2- tdpdfor) as the reference level
msm$ForagingType <- factor(msm$ForagingType, levels = c(2, 1, 3, 4))

tmod <- MASS::glm.nb(tot_affil ~ offset(log(avail_assoc)) +
                                 matkins + patkins +
                                 distkins +
                                 ForagingType + provisioned + adult_age +
                                 sightings + adult_yrs,
                    data = msm)

umod <- MASS::glm.nb(unrel_affil ~ offset(log(avail_assoc)) +
                                   matkins + patkins +
                        distkins +
                                   ForagingType + provisioned + adult_age +
                                  sightings + adult_yrs,
                     data = msm)

summary(tmod)

t_est <- as.data.frame(cbind(Estimate = coef(tmod), confint(tmod)))
t_est$pv <- ifelse(sign(t_est[, 2]) == sign(t_est)[, 3], "*", "")

t_est

summary(umod)

u_est <- as.data.frame(cbind(Estimate = coef(umod), confint(umod)))
u_est$pv <- ifelse(sign(u_est[, 2]) == sign(u_est)[, 3], "*", "")

u_est

########################################
## Centrality Data
########################################

rnms <- real_network_metrics
names(rnms)[-1] <- paste0("raw_", names(rnms)[-1])
# get median of all simulated metrics
arm <- do.call("rbind", central_rand)
med_random_metrics <- aggregate(. ~ ego, data = arm, median)
med_random_metrics <- med_random_metrics[match(rnms$ego, med_random_metrics$ego), ]
names(med_random_metrics) <- c(ego, paste0("sim_", names(med_random_metrics)[-1]))
all(rnms$ego == med_random_metrics$ego)

# create residuals
std_s_metrics <- as.data.frame(as.matrix(rnms[, 2:5]) - as.matrix(med_random_metrics[, 2:5]))
names(std_s_metrics) <- gsub("raw", "resid", names(std_s_metrics))

sm <- cbind(rnms, std_s_metrics, med_random_metrics)

# add in network metrics
msm <- merge(msm, sm, by.x = "Dolphin.ID", by.y = "ego")

# model centrality
emod <- lm(scale(resid_eigen) ~ matkins + patkins +
            distkins +
            ForagingType + provisioned +
            adult_age + sightings + adult_yrs, data = msm)

summary(emod)
car::vif(emod)

e_est <- as.data.frame(cbind(Estimate = coef(emod), confint(emod)))
e_est$pv <- ifelse(sign(e_est[, 2]) == sign(e_est)[, 3], "*", "")

e_est

cmod <- lm(scale(resid_closeness) ~ matkins + patkins +
            distkins +
            ForagingType + provisioned +
            adult_age + sightings + adult_yrs, data = msm)

summary(cmod)
car::vif(cmod)

c_est <- as.data.frame(cbind(Estimate = coef(cmod), confint(cmod)))
c_est$pv <- ifelse(sign(c_est[, 2]) == sign(c_est)[, 3], "*", "")

c_est

# save(tmod, t_est, umod, u_est, emod, e_est, cmod, c_est,
#      msm, msm,
#      file = "IntermediateData/model_results.RData")

# write out table 2 info

mod_results <- t_est

rownames(mod_results) <- c("Intercept",
                       "Close Matrilineal Kin",
                       "Close Non-Matrilineal Kin",
                       "Distant Kin",
                       "Foraging: Mixed",
                       "Foraging: Sponging",
                       "Foraging: Seagrass",
                       "Provisioning",
                       "Age (yrs)",
                       "log(Number Observations)",
                       "Years Observed")

results <- cbind(exp(mod_results[, 1]), mod_results)
results <- cbind(apply(results[, 1:4], 2, function(x) round(x, 3)), results[, 5])
colnames(results) <- c("Odds Ratio", "Beta", "Lower", "Upper", "Significant")
finalresults <- data.frame(results)
finalresults$Beta <- paste0(finalresults$Beta, " (",
                                  finalresults$Lower, " - ",
                                  finalresults$Upper, ")",
                            finalresults$Significant)

mod_results <- u_est

results <- cbind(exp(mod_results[, 1]), mod_results)
results <- cbind(apply(results[, 1:4], 2, function(x) round(x, 3)), results[, 5])
colnames(results) <- c("Odds Ratio Unaffiliated", "Beta_Unaffiliated", "Lower", "Upper", "Significant")
results <- data.frame(results)
results$Beta_Unaffiliated <- paste0(results$Beta, " (",
                              results$Lower, " - ",
                              results$Upper, ")",
                              results$Significant)

results <- results[, 1:2]


finalresults <- cbind(finalresults, results)
finalresults <- subset(finalresults, select = -c(Lower, Upper, Significant))

# write.csv(finalresults, file="TablesFigures/Table2_affiliations.csv")

# write out table 3 info
mod_results <- e_est

rownames(mod_results) <- c("Intercept",
                           "Close Matrilineal Kin",
                           "Close Non-Matrilineal Kin",
                           "Distant Kin",
                           "Foraging: Mixed",
                           "Foraging: Sponging",
                           "Foraging: Seagrass",
                           "Provisioning",
                           "Age (yrs)",
                           "log(Number Observations)",
                           "Years Observed")

results <- apply(mod_results[-ncol(mod_results)], 2, function(x) round(x, 4))
colnames(results) <- c("Estimate", "Lower", "Upper")
finalresults <- data.frame(results, Significant = mod_results$pv)
finalresults$Estimate <- paste0(finalresults$Estimate, " (",
                                  finalresults$Lower, " - ",
                                  finalresults$Upper, ")")

mod_results <- c_est
results <- apply(mod_results[-ncol(mod_results)], 2, function(x) round(x, 4))
colnames(results) <- c("Estimate_Closeness", "Lower", "Upper")
results <- data.frame(results, Significant = mod_results$pv)
results$Estimate_Closeness <- paste0(results$Estimate_Closeness, " (",
                             results$Lower, " - ",
                             results$Upper, ")")

results <- results[, -(2:3)]

finalresults <- cbind(finalresults, results)
finalresults <- subset(finalresults, select = -c(Lower, Upper))

# write.csv(finalresults, file="TablesFigures/Table3_centrality.csv")
