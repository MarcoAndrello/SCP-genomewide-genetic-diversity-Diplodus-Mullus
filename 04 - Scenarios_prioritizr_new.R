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

for (i.species in 1 : 2) {
    
    name_species = c("Diplodus_sargus", "Mullus_surmuletus")[i.species]
    num_axes <- c(17, 26)[i.species]
    
    # Create an object containing the coordinates of the occupied PU
    # in all genetic PCA axes: need this for finding clusters
    pus_g_values <- cbind(g_rast_values[[i.species]], st_drop_geometry(pus))
    names(pus_g_values)[1:ncol(g_rast_values[[i.species]])] <- paste0("pca",sprintf("%02d",1:ncol(g_rast_values[[i.species]])))
    species_coord <- filter(pus_g_values, !!as.name(name_species) == 1)
    species_coord <- species_coord[1:num_axes]
    rm(pus_g_values)
    
    # Run kmeans with several number of clusters
    set.seed(123)
    wss <- vector()
    for (i.clust in 1 : 20) {
        clusters <- kmeans(species_coord,centers=i.clust)
        wss[i.clust] <- clusters$tot.withinss
    }
    
    # Plot k-means
    png(paste0("Number_of_clusters_",name_species,".png"),width=10,height=10,units="cm",res=300)
    plot(c(1:20),wss,type="l",xlab="Number of clusters",ylab="Within sum-of-squares",
         main = name_species)
    dev.off()
}

# For Diplodus: it suggests 5 (actually 7!!) clusters. We use K = 2 (Boulanger et al 2022) and K = 5
# For Mullus: it suggests 7 clusters. We use K = 3 (Boulanger et al 2022) and K = 7

# mi fermo qui

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
        for (i.axis in 1 : 17) {
            st_matrix <- split.taxon(x = g_rast_values_Diplodus[,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,"Diplodus_sargus"),
                                     name_feat = paste0("Diplodus_a",sprintf("%02d",i.axis)))
            pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
        }
        # Mullus surmuletus
        for (i.axis in 1 : 26) {
            st_matrix <- split.taxon(x = g_rast_values_Mullus[,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,"Mullus_surmuletus"),
                                     name_feat = paste0("Mullus_a",sprintf("%02d",i.axis)))
            pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
        }
    
        
        
        

        # Define classes
        pus_split.taxon <- pus
        for (i.axis in 1 : num_axes) {
            # Split-taxon function
            st_matrix <- split.taxon(x = g_rast_values[,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,name_species),
                                     name_feat = paste0(name_species_feat,sprintf("%02d",i.axis)))
            pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
        }
        ## Solve problem
        features <- names(pus_split.taxon)[grep(name_species_feat,names(pus_split.taxon))]
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

################################################################################
# Combining PCA axes to define conservation features
################################################################################
prob_multi <- res_multi <- list()
for (i.num_classes in 1 : 2){
    num_classes <- v.num_classes.multi[i.num_classes]
    st_matrix <- split.taxon.multi(x = g_rast_values[,1:num_axes],
                                   num_classes = num_classes,
                                   rij_taxon = pull(pus,name_species),
                                   name_feat = paste0(species,"_m",num_axes))
    pus_split.taxon <- cbind(pus, st_matrix)
    ## Solve problem
    features <- names(pus_split.taxon)[grep(paste0(species,"_m",num_axes),names(pus_split.taxon))]
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
      paste0("multi_kmeans_",v.num_classes.multi))

save(problems,
     results,
     file=paste0(getwd(),"/Results/Results_prioritizr_",species,".RData"))
