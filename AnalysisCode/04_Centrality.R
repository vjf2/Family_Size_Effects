#######################################
#### Calculate Adult Network Centrality
#######################################

library(igraph)
library(SocGen)
library(foreach)
library(doParallel)

indiv_covars <- read.csv("SharedData/individual_covariates.csv")

focal_adults <- indiv_covars[which(indiv_covars$num_sightings_adult >= 20 &
                                     indiv_covars$relocations >= 35 &
                                     indiv_covars$genotyped == "Y" &
                                     indiv_covars$prop_genotyped >= 0.70 &
                                     !is.na(indiv_covars$adult_entry) &
                                     !is.na(indiv_covars$provisioned)), ]

# sort these females into groups by available years to do calculations of the network while they are in their adult period

sightings <- read.csv("SharedData/sightings.csv")

dates <- sort(unique(sightings$Date))

focal_adults$first_obs_date <- as.Date(sapply(1:length(focal_adults$Dolphin.ID),
                                            function(x) {dates[dates >= focal_adults$adult_entry[x]][1]}),
                                     origin = "1970-01-01")

focal_adults$last_obs_date <- as.Date(sapply(1:length(focal_adults$Dolphin.ID), function(x) {
  tail(dates[dates <= focal_adults$depart[x]], 1)
}), origin = "1970-01-01")

# date range key table

date_range <- unique(focal_adults[, c("first_obs_date", "last_obs_date")])

# add date diff
date_range$datediff <- as.numeric((date_range$last_obs_date - date_range$first_obs_date) / 365.25)

date_range$date_key <- 1:nrow(date_range)

focal_adults <- merge(focal_adults, date_range, sort = FALSE)

# set up empty dataframe

nf <- nrow(focal_adults)

real_network_metrics <- data.frame(ego = focal_adults$Dolphin.ID,
                                 degree = numeric(nf),
                                 strength = numeric(nf),
                                 eigen = numeric(nf),
                                 closeness = numeric(nf),
                                 network_size = numeric(nf),
                                 gc_size = numeric(nf))

for (i in 1:nrow(date_range)) {

  centrality_surveys <- sightings[sightings$Date >= date_range$first_obs_date[i] &
                                  sightings$Date <= date_range$last_obs_date[i], ]

  central_network <- half_weight(sightings = centrality_surveys,
                               group_variable = "Observation.ID",
                               dates = "Date",
                               IDs = "Dolphin.ID")

  g <- graph.adjacency(central_network, mode = "undirected", weighted = TRUE, diag = FALSE)

  network_size <- vcount(g)

  graphs <- decompose.graph(g)
  largest <- which.max(sapply(graphs, vcount))

  g1 <- graphs[[largest]]
  gc_size <- vcount(g1)

  if (i == 3) {print(distance_table(g)); print(mean_distance(g))} # whole network

  egos <- focal_adults$Dolphin.ID[focal_adults$date_key == i]

  # mini_loop
  for (j in egos) {
    # calc network metrics here
    real_network_metrics[real_network_metrics$ego == j, "degree"] <- degree(g, j)
    real_network_metrics[real_network_metrics$ego == j, "strength"] <- strength(g, j)
    real_network_metrics[real_network_metrics$ego == j, "eigen"] <- eigen_centrality(g)$vector[j]
    real_network_metrics[real_network_metrics$ego == j, "closeness"] <- closeness(g,
                                                                                v = j,
                                                                                weights = 1 / E(g)$weight,
                                                                                normalized = TRUE)
    # report network size and giant component size
    real_network_metrics[real_network_metrics$ego == j, "network_size"] <- network_size
    real_network_metrics[real_network_metrics$ego == j, "gc_size"] <- gc_size

    # note: closeness centrality is not well-defined for disconnected graphs

  }
}

# repeat for random
load("IntermediateData/kfinal1000.RData")

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {library(SocGen);library(igraph)})
clusterExport(cl, c("kfinal", "date_range", "focal_adults", "nf"))
registerDoParallel(cl)

starttime <- Sys.time()

central_rand <- foreach(n = 1:length(kfinal), .errorhandling = 'pass') %dopar% {

  kfinal1 <- kfinal[[n]]
  kfinal1$Date <- substring(kfinal1$group, 1, 10)

  network_metrics <- data.frame(ego = focal_adults$Dolphin.ID,
                              degree = numeric(nf),
                              strength = numeric(nf),
                              eigen = numeric(nf),
                              closeness = numeric(nf))

  for (i in 1:nrow(date_range)) {

    centrality_surveys <- kfinal1[kfinal1$Date >= date_range$first_obs_date[i] &
                                  kfinal1$Date <= date_range$last_obs_date[i], ]

    central_network <- half_weight(sightings = centrality_surveys,
                                 group_variable = "group",
                                 dates = "Date",
                                 IDs = "id")

    g <- graph.adjacency(central_network, mode = "undirected", weighted = TRUE, diag = FALSE)

    egos <- focal_adults$Dolphin.ID[focal_adults$date_key == i]

    # mini_loop
    for (j in egos) {
      # calc network metrics here
      network_metrics[network_metrics$ego == j, "degree"] <- degree(g, j)
      network_metrics[network_metrics$ego == j, "strength"] <- strength(g, j)
      network_metrics[network_metrics$ego == j, "eigen"] <- eigen_centrality(g)$vector[j]
      network_metrics[network_metrics$ego == j, "closeness"] <- closeness(g,
                                                                        v = j,
                                                                        weights = 1 / E(g)$weight,
                                                                        normalized = TRUE)

      # note: closeness centrality is not well-defined for disconnected graphs
    }

    }
  # cat(paste0(" network complete for number ", n, "\n")) #if printing to file

  network_metrics
}

endtime <- Sys.time()

stopCluster(cl)

endtime - starttime # run time ~1 hour

beepr::beep(2)

# save relevant objects

focal_adults <- focal_adults[, names(indiv_covars)]

save(real_network_metrics, central_rand, focal_adults, file = "IntermediateData/centrality_networks.RData")
