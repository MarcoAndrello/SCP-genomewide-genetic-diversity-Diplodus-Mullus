# Scenarios prioritizr
rm(list=ls())

species <- "Mullus"

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)
library(parallel)
library(rgeoda)

if (species == "Diplodus") {
    name_species = "Diplodus_sargus"
    genetic_raster = "/Results May_2023/Diplodus_allAxes_8068.grd"
    name_species_feat = "Diplodus_a"
    num_axes = 17
    v.num_classes.multi = c(2,5)
}
if (species == "Mullus") {
    name_species = "Mullus_surmuletus"
    genetic_raster = "/Results May_2023/Mullus_allAxes_2753.grd"
    name_species_feat = "Mullus_a"
    num_axes = 26
    v.num_classes.multi = c(3,7)
}

# Number of assigned threads
threads <- 1L #as.integer(max(1, detectCores() - 2))

# Source functions
source(paste0(getwd(),"/functions/split.taxon.R"))

# Load data
## Planning units
load("Planning_units.RData")
## Genetic data
g_rast <- rast(paste0(getwd(),genetic_raster))

# Extract genetic values from raster for the centroid of each PU
g_rast_values <- terra::extract(g_rast,pus_centroid)

# Create an object containing the coordinates of the occupied PU
# in all genetic PCA axes: need this for finding clusters
pus_g_values <- cbind(g_rast_values, st_drop_geometry(pus))
names(pus_g_values)[1:ncol(g_rast_values)] <- paste0("pca",sprintf("%02d",1:ncol(g_rast_values)))
species_coord <- filter(pus_g_values, !!as.name(name_species) == 1)
species_coord <- species_coord[1:num_axes]
rm(pus_g_values)
# Run kmeans with several number of clusters
wss <- vector()
for (i.clust in 1 : 20) {
    clusters <- kmeans(species_coord,centers=i.clust)
    wss[i.clust] <- clusters$tot.withinss
}
png(paste0("Number_of_clusters_",species,".png"),width=10,height=10,units="cm",res=300)
plot(c(1:20),wss,type="l",xlab="Number of clusters",ylab="Within sum-of-squares",
     main = name_species)
dev.off()
# For Diplodus: it suggests 5 clusters. We use K = 2 (Boulanger et al 2022) and K = 5
# For Mullus: it suggests 7 clusters. We use K = 3 (Boulanger et al 2022) and K = 7


################################################################################
# Using single PCA axes as conservation features
################################################################################

v.class_method <- c("quantile","equal","natural","stddev")
v.num_classes <- c(3,6,12)
prob_single <- res_single <- list()
i.scen <- 1
for (i.class_method in 1 : 4){
    class_method <- v.class_method[i.class_method]
    for (i.num_classes in 1 : 3){
        num_classes <- v.num_classes[i.num_classes]
        if (class_method == "stddev" & num_classes != 6) next # only 6 classes for stddev
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
    p <- problem(pus_split.taxon,
                 features = features,
                 cost_column = "cost") %>%
        add_min_set_objective() %>%
        add_relative_targets(0.15) %>%
        add_binary_decisions() %>%
        add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
        add_gurobi_solver(gap=0.02, threads=threads)
    prob_multi[[i.num_classes]] <- p
    res_multi[[i.num_classes]] <- solve(p, run_checks=F)
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
     file=paste0("Results_prioritizr_",species,".RData"))
