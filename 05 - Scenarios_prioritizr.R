# Scenarios prioritizr
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)
library(parallel)
library(rgeoda)

# Number of assigned threads
threads <- 1L #as.integer(max(1, detectCores() - 2))

# Source functions
source(paste0(getwd(),"/functions/split.taxon.R"))

# Load data
## Planning units
load("Planning_units.RData")

# Open genetic raster and Extract values for the centroid of each PU
g_rast_values <- list()
g_rast <- rast(paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))
g_rast_values[[1]] <- terra::extract(g_rast,pus_centroid)
g_rast <- rast(paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))
g_rast_values[[2]] <- terra::extract(g_rast,pus_centroid)
rm(g_rast)


################################################################################
# Using single PCA axes as conservation features
################################################################################

v.class_method <- c("quantile","equal","natural","stddev")
v.num_classes <- c(3,6,12)
prob_single <- res_single <- list()
i.scen <- 1
i.class_method <- i.num_classes <- 1
for (i.class_method in 1 : 4){
    class_method <- v.class_method[i.class_method]
    for (i.num_classes in 1 : 3){
        num_classes <- v.num_classes[i.num_classes]
        if (class_method == "stddev" & num_classes != 6) next # only 6 classes for stddev
        
        # Diplodus sargus
        pus %>% 
            select(Diplodus_sargus, Mullus_surmuletus, cost, status) -> 
            pus_split.taxon
        for (i.axis in 1 : 14) {
            st_matrix <- split.taxon(x = g_rast_values[[1]][,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,"Diplodus_sargus"),
                                     name_feat = paste0("Diplodus_a",sprintf("%02d",i.axis)))
            pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
        }
        # Mullus surmuletus
        for (i.axis in 1 : 21) {
            st_matrix <- split.taxon(x = g_rast_values[[2]][,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,"Mullus_surmuletus"),
                                     name_feat = paste0("Mullus_a",sprintf("%02d",i.axis)))
            pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
        }
        # Feature names
        features <- names(pus_split.taxon)[c(grep("Diplodus_a",names(pus_split.taxon)),
                                             grep("Mullus_a",names(pus_split.taxon)))]

        ## Solve problem
        prob_single[[i.scen]] <- problem(pus_split.taxon,
                     features = features,
                     cost_column = "cost") %>%
            add_min_set_objective() %>%
            add_relative_targets(0.15) %>%
            add_binary_decisions() %>%
            add_locked_in_constraints(which(pus_split.taxon$status==1)) %>%
            add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
            add_gurobi_solver(gap=0.02, threads=threads)
        res_single[[i.scen]] <- solve(prob_single[[i.scen]], run_checks=F)
        i.scen <- i.scen + 1
    }
}
rm(st_matrix,class_method,num_classes,i.class_method,i.num_classes,i.scen,i.axis,v.class_method)


#### QUI DA MODIFICARE DANDO A SPLIT MULTI GIA LA CLUSTERIZZAZIONE FATTA
################################################################################
# Combining PCA axes to define conservation features
################################################################################
prob_multi <- res_multi <- list()
for (i.num_classes in 1 : 2){
    pus %>% 
        select(Diplodus_sargus, Mullus_surmuletus, cost, status, cost) -> 
        pus_split.taxon
    # Diplodus
    num_classes <- c(2,7)[i.num_classes] #c(2,5)[i.num_classes]
    st_matrix <- split.taxon.multi(x = g_rast_values[[1]][,1:17],
                                   num_classes = num_classes,
                                   rij_taxon = pull(pus,"Diplodus_sargus"),
                                   name_feat = paste0("Diplodus_m14"))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
    # Mullus
    num_classes <- c(3,7)[i.num_classes]
    st_matrix <- split.taxon.multi(x = g_rast_values[[2]][,1:26],
                                   num_classes = num_classes,
                                   rij_taxon = pull(pus,"Mullus_surmuletus"),
                                   name_feat = paste0("Mullus_m21"))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
    # Feature names
    features <- names(pus_split.taxon)[c(grep("Diplodus_m",names(pus_split.taxon)),
                                         grep("Mullus_m",names(pus_split.taxon)))]
    ## Solve problem
    prob_multi[[i.num_classes]] <- problem(pus_split.taxon,
                                           features = features,
                                           cost_column = "cost") %>%
        add_min_set_objective() %>%
        add_relative_targets(0.15) %>%
        add_binary_decisions() %>%
        add_locked_in_constraints(which(pus_split.taxon$status==1)) %>%
        add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
        add_gurobi_solver(gap=0.02, threads=threads)
    res_multi[[i.num_classes]] <- solve(prob_multi[[i.num_classes]], run_checks=F)
}

problems <- c(prob_single,prob_multi)
results <- c(res_single,res_multi)
names(problems) <- names(results) <- 
    c("single_quantile_3",
      "single_quantile_6",
      "single_quantile_12",
      "single_equal_3",
      "single_equal_6",
      "single_equal_12",
      "single_natural_3",
      "single_natural_6",
      "single_natural_12",
      "single_sd_6",
      "multi_kmeans_low",
      "multi_kmeans_high")

save(problems,
     results,
     file=paste0(getwd(),"/Results/Results_prioritizr.RData"))









# For Diplodus: it suggests 5 (actually 7!!) clusters. We use K = 2 (Boulanger et al 2022) and K = 5
# For Mullus: it suggests 7 clusters. We use K = 3 (Boulanger et al 2022) and K = 7

# Run kmeans with several number of clusters
set.seed(123)
wss <- vector()
for (i.clust in 1 : 20) {
    # clusters <- kmeans(species_coord,centers=i.clust)
    pam[[i.clust]] <- pam(species_coord,i.clust)
    # wss[i.clust] <- clusters$tot.withinss
}

# Plot k-means
png(paste0("Number_of_clusters_",name_species,".png"),width=10,height=10,units="cm",res=300)
plot(c(1:20),wss,type="l",xlab="Number of clusters",ylab="Within sum-of-squares",
     main = name_species)
dev.off()