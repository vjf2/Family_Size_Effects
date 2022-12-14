#######################################
#### Calculate Association Coefficients
#######################################

library(SocGen)
library(igraph)

load("IntermediateData/kfinal1000.RData")

indiv_covars <- read.csv("SharedData/individual_covariates.csv")
indiv_covars$entry <- as.Date(indiv_covars$entry)
indiv_covars$depart <- as.Date(indiv_covars$depart)

sightings <- read.csv("SharedData/sightings.csv")

dates <- sort(unique(sightings$Date))

# select pairs to from which to determine affiliations
# affiliated pairs are calculated between all individuals with at least 2 years postweaning, 35 relocations and 20 sightings in which group associations were recorded and who have been genotyped

affil_females <- indiv_covars[which(indiv_covars$total_yrs >= 2 &
                                    indiv_covars$genotyped == "Y" &
                                    indiv_covars$relocations >= 35 &
                                    indiv_covars$num_sightings >= 20), ]

# mask the data so that association rates are only estimated in the timeframe where both members are alive
affil_mask <- schedulize(affil_females,
                         id = "Dolphin.ID",
                         start = "entry",
                         end = "depart",
                         dates = dates,
                         format = "mask")

affil_sightings <- sightings[sightings$Dolphin.ID %in% affil_females$Dolphin.ID, ]

masked_network <- half_weight(sightings = affil_sightings,
                              group_variable = "Observation.ID",
                              dates = "Date",
                              IDs = "Dolphin.ID",
                              mask = affil_mask)

# repeat for the results of the random model

library(foreach)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(SocGen)})
clusterExport(cl, c("kfinal", "affil_mask", "affil_females"))
registerDoParallel(cl)

starttime <- Sys.time()

ai_affil_rand <- foreach(n = 1:length(kfinal), .errorhandling = 'pass') %dopar% {

  kfinal1 <- kfinal[[n]]
  kfinal1$Date <- substring(kfinal1$group, 1, 10)

  affil_kfinal1 <- kfinal1[kfinal1$id %in% affil_females$Dolphin.ID, ]

  masked_network <- half_weight(sightings = affil_kfinal1,
                              group_variable = "group",
                              dates = "Date",
                              IDs = "id",
                              mask = affil_mask)

  cat(paste0(" network complete for number ", n, "\n"))

  masked_network
}

endtime <- Sys.time()

stopCluster(cl)

endtime - starttime # expected run time ~3 min


save(masked_network, ai_affil_rand, file = "IntermediateData/affiliation_networks.RData")
