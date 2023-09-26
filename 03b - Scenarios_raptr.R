# Scenarios raptr
rm(list=ls())

species <- "Mullus"

library(tidyverse)
library(sf)
library(terra)
library(raptr)
library(parallel)
library(rgeoda)
library(DescTools)

if (species == "Diplodus") {
    name_species = "Diplodus_sargus"
    genetic_raster = "/Results May_2023/Diplodus_allAxes_8068.grd"
    name_species_feat = "Diplodus_a"
    num_axes = 17
}
if (species == "Mullus") {
    name_species = "Mullus_surmuletus"
    genetic_raster = "/Results May_2023/Mullus_allAxes_2753.grd"
    name_species_feat = "Mullus_a"
    num_axes = 26
}

# Number of assigned threads
threads <- 8L #as.integer(max(1, detectCores() - 2))

# Source functions
source(paste0(getwd(),"/functions/split.taxon.R"))

# Load data
## Planning units
load("Planning_units.RData")
pus$status <- as.integer(pus$status) # status should be integer
## Genetic data
g_rast <- rast(paste0(getwd(),genetic_raster))

# Extract genetic values from raster for the centroid of each PU
g_rast_values <- terra::extract(g_rast,pus_centroid)

# Create an object containing the coordinates of the occupied PU
# in all genetic spaces (i.e. for all genetic PCA axes)(need this for planning.unit.points)
pus_g_values <- cbind(g_rast_values, st_drop_geometry(pus))
names(pus_g_values)[1:ncol(g_rast_values)] <- paste0("pca",sprintf("%02d",1:ncol(g_rast_values)))
species_coord <- filter(pus_g_values, !!as.name(name_species) == 1) 
species_coord <- species_coord[1:ncol(g_rast_values)]
rm(pus_g_values)

# Create species probabilities
pu.species.probabilities <- data.frame(species=1L,
                                       pu=c(1:nrow(pus)),
                                       value=pus[[name_species]])
pu.species.probabilities <- filter(pu.species.probabilities, value==1)

# Coonvert to PolySet and calculate boundary
pus %>%  convert2PolySet() -> polygons
boundary <- raptr::calcBoundaryData(polygons)


################################################################################
# Using single PCA axes as genetic spaces
################################################################################
v.class_method <- c("quantile","equal","natural","stddev")
v.num_classes <- c(3,6,12)
attribute.spaces <- list() # contains elements = axes of the PCA
prob_single <- res_single <- list()
i.scen <- 1
for (i.class_method in 1 : 4){
    class_method <- v.class_method[i.class_method]
    for (i.num_classes in 1 : 3){
        num_classes <- v.num_classes[i.num_classes]
        if (class_method == "stddev" & num_classes != 6) next # only 6 classes for stddev
        # Define classes
        for (i.axis in 1 : num_axes) {
            res_split.taxon <- split.taxon(x = g_rast_values[,i.axis],
                                     num_classes = num_classes,
                                     class_method = class_method,
                                     rij_taxon = pull(pus,name_species),
                                     name_feat = paste0(name_species_feat,sprintf("%02d",i.axis)),
                                     return_class_midpoint = T)
            coords <- as.matrix(res_split.taxon[[1]])
            st_matrix <- res_split.taxon[[2]]
            weights <- colSums(st_matrix)/sum(st_matrix) # Weighting DP by the proportion of PU in the class they represent
            # if a dp has 0 weight, remove it
            id0 <- which(weights==0)
            if (length(id0) > 0) {
                coords <- matrix(coords[-id0,],ncol=1)
                weights <- weights[-id0]
                cat("remove",length(id0),"DP from",class_method,num_classes,"\n")
            }
            # Demand points
            dp <- DemandPoints(coords, weights)
            # Attribute space
            genetic_spaces_Axis <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=as.matrix(species_coord[,i.axis]),
                                                                                                 ids=pu.species.probabilities$pu), # only occupied PUs
                                                       demand.points = dp,
                                                       species=1L))
            # Attribute spaces for this axis
            attribute.spaces[[i.axis]] <-
                AttributeSpaces(
                    spaces = genetic_spaces_Axis,
                    name = paste0("genetic_",i.axis)
                )
            
        }
        # Create targets (Geographic targets are set to 0)
        targets <- data.frame(species = 1L,
                              target = as.integer(0:num_axes),
                              proportion = c(0.15,rep(0.75,num_axes))
        )
        # Create rap_data object
        rap_data <- RapData(pu = st_drop_geometry(pus),
                            species = data.frame(name=name_species),
                            targets = targets,
                            pu.species.probabilities = pu.species.probabilities,
                            attribute.spaces = attribute.spaces,
                            boundary = boundary,
                            polygons = polygons)
        # Define problem and solve it
        ro <- RapUnreliableOpts(BLM=0)
        prob_single[[i.scen]] <- RapUnsolved(ro, rap_data)
        res_single[[i.scen]] <- solve(prob_single[[i.scen]],
                                      Threads = threads, verbose=T, NumericFocus= 3L, MIPGap=0.02, NumberSolutions=100L)
        
        i.scen <- i.scen + 1
    }
}
save(prob_single,res_single,file=paste0("Results_raptr_Mullus_PARTIAL_",i.scen,".RData"))


################################################################################
# Combining PCA axes to define genetic spaces
################################################################################
# # Using kmeans
# v.num_classes <- c(3,6,12)
# attribute.spaces <- list()
# prob_multi <- res_multi <- list()
# i.num_classes <- 1
# for (i.num_classes in 1 : 3){
#     num_classes <- v.num_classes[i.num_classes]
#     res_split.taxon <- split.taxon.multi(x = g_rast_values[,1:num_axes],
#                                    num_classes = num_classes,
#                                    rij_taxon = pull(pus,name_species),
#                                    name_feat = paste0(species,"_m",num_axes),
#                                    return_class_midpoint = T)
#     coords <- as.matrix(res_split.taxon[[1]]) # Centroids of each clusters in the num_axes dimensions
#     st_matrix <- res_split.taxon[[2]]
#     weights <- colSums(st_matrix)/sum(st_matrix) # Weighting DP by the proportion of PU in the class they represent
#     id0 <- which(weights==0)
#     if (length(id0) > 0) {
#         coords <- matrix(coords[-id0,],ncol=1)
#         weights <- weights[-id0]
#         cat("remove",length(id0),"DP from",num_classes,"\n")
#     }
#     dp <- DemandPoints(coords, weights)
#     genetic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=as.matrix(species_coord[,1:num_axes]),
#                                                                                     ids=pu.species.probabilities$pu), # only occupied PUs
#                                           demand.points = dp,
#                                           species=1L))
#     attribute.spaces[[1]] <-
#         AttributeSpaces(
#             spaces = genetic_spaces,
#             name = "genetic"
#         )
#     # Create targets
#     targets <- data.frame(species = 1L,
#                           target = as.integer(0:1),
#                           proportion = c(0.15,0.75)
#     )
#     
#     rap_data <- RapData(pu = st_drop_geometry(pus),
#                         species = data.frame(name=name_species),
#                         targets = targets,
#                         pu.species.probabilities = pu.species.probabilities,
#                         attribute.spaces = attribute.spaces,
#                         boundary = boundary,
#                         polygons = polygons)
#     # Define problem and solve it
#     ro <- RapUnreliableOpts(BLM=0)
#     prob_multi[[i.num_classes]] <- RapUnsolved(ro, rap_data)
#     print(maximum.targets(prob_multi[[i.num_classes]]))
#     res_multi[[i.num_classes]] <- solve(prob_multi[[i.num_classes]],
#                                          Threads = threads, verbose=T, NumericFocus= 3L, MIPGap=0.02, NumberSolutions=100L)
# }
##

# Using hypervolume
attribute.spaces <- list()
prob_multi_hypervolume <- res_multi_hypervolume <- list()
v.num_dp <- c(20,40,80)
i.num_dp <- 1
for (i.num_dp in 1 : 3) {
    g_rast_values[,c(1:num_axes)] %>%
        filter(pus$Diplodus_sargus==1) %>%
        as.matrix() %>%
        unique() -> # remove duplicate rows (there are many (bcs of the spatial interpolation method), and speeds up the calculation)
        points              #
    dp <- make.DemandPoints(points, quantile = 0.95,
                            kernel.method = "hypervolume", samples.per.point=1000, n=as.integer(v.num_dp[i.num_dp]))
    # Attribute space       # In this list there will be AttributeSpace for each species
    genetic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=as.matrix(species_coord[,1:num_axes]),
                                                                                    ids=pu.species.probabilities$pu), # only occupied PUs
                                          demand.points = dp,
                                          species=1L))
    # Attribute spacews for this axis (both species)
    attribute.spaces[[1]] <-
        AttributeSpaces(
            spaces = genetic_spaces,
            name = "genetic"
        )
    # Create targets
    ## Geographic targets to 0
    targets <- data.frame(species = 1L,
                          target = as.integer(0:1),
                          proportion = c(0.15,0.75)
    )
    rap_data <- RapData(pu = st_drop_geometry(pus),
                        species = data.frame(name=species),
                        targets = targets,
                        pu.species.probabilities = pu.species.probabilities,
                        attribute.spaces = attribute.spaces,
                        boundary = boundary,
                        polygons = polygons)
    # Define problem and solve it
    ro <- RapUnreliableOpts(BLM=0)
    prob_multi_hypervolume[[i.num_dp]] <- RapUnsolved(ro, rap_data)
    print(maximum.targets(prob_multi_hypervolume[[i.num_dp]]))
    res_multi_hypervolume[[i.num_dp]] <- solve(prob_multi_hypervolume[[i.num_dp]],
                                               Threads = threads, verbose=T, NumericFocus= 3L, MIPGap=0.02, NumberSolutions=100L)
}
save(prob_multi_hypervolume,res_multi_hypervolume,file="Results_raptr_Mullus_PARTIAL_HYPERVOLUME.RData")

###

problems <- c(prob_single, prob_multi_hypervolume)
results <- c(res_single, res_multi_hypervolume)
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
     "multi_hypervolume_20",
     "multi_hypervolume_40",
     "multi_hypervolume_80")
save(problems,
     results,
     file=paste0("Results_raptr_",species,".RData")
)
