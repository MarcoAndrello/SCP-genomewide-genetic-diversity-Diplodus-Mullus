# Evaluation and extension of the current set of MR
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(raptr)
library(prioritizr)
library(parallel)
library(tictoc)
library(tmap)
library(rnaturalearth)
library(RColorBrewer)

# Number of assigned threads
threads <- 8L #as.integer(max(1, detectCores() - 2))

# Load data
## Planning units
load("Planning_units.RData")

################################################################################
# Definition of raptr problem for both species
################################################################################

# Create pu.species.probabilities
pu.species.probabilities <- rbind(data.frame(species=1L,
                                             pu=c(1:nrow(pus)),
                                             value=pus[["Diplodus_sargus"]]),
                                  data.frame(species=2L,
                                             pu=c(1:nrow(pus)),
                                             value=pus[["Mullus_surmuletus"]]))
pu.species.probabilities <- filter(pu.species.probabilities, value==1)

# Create targets
## Geographic targets to 0
targets <- rbind(data.frame(species = 1L,
                            target = as.integer(0:1),
                            proportion = c(0.15,0.95)),
                 data.frame(species = 2L,
                            target = as.integer(0:1),
                            proportion = c(0.15,0.95))
                 )
                            
genetic_spaces <- list()
for (i.species in 1 : 2) {
    
    species <- c("Diplodus","Mullus")[i.species]
    
    if (species == "Diplodus") {
        name_species = "Diplodus_sargus"
        genetic_raster = "/Results/Diplodus_allAxes_8068.grd"
        name_species_feat = "Diplodus_a"
        num_axes = 17
    }
    if (species == "Mullus") {
        name_species = "Mullus_surmuletus"
        genetic_raster = "/Results/Mullus_allAxes_2753.grd"
        name_species_feat = "Mullus_a"
        num_axes = 26
    }
    ## Genetic data
    g_rast <- rast(paste0(getwd(),genetic_raster))
    
    # Extract genetic values from raster for the centroid of each PU
    g_rast_values <- terra::extract(g_rast,pus_centroid)
    
    # Create an object containing the coordinates of the occupied PU
    # in all genetic spaces (i.e. for all genetic PCA axes)(need this for planning.unit.points)
    pus_g_values <- cbind(g_rast_values, st_drop_geometry(pus))
    names(pus_g_values)[1:ncol(g_rast_values)] <- paste0("pca",sprintf("%02d",1:ncol(g_rast_values)))
    species_coord <- filter(pus_g_values, !!as.name(name_species) == 1) 
    species_coord <- species_coord[1:num_axes]
    rm(pus_g_values)
    
    # Many PUs MIGHT have the same multidimensional coordinates in the PCA space
    # (this happens only with the nearest neighbor interpolation method)
    # Count number of PUs per multidimensional coordinate
    species_coord %>% group_by_all() %>% count -> coords
    # The number of PUs per coordinate will be the weight of that coordinate
    coords[,num_axes+1] %>% pull(n) -> weights
    weights <- weights/sum(weights)
    # Remove last column (n, the counts)
    coords <- coords[,1:num_axes]
    # Make demand points
    dp <- DemandPoints(as.matrix(coords), weights)
    
    # Attribute space       # In this list there will be AttributeSpace for each species
    genetic_spaces[[i.species]] <- AttributeSpace(planning.unit.points = PlanningUnitPoints(coords = as.matrix(species_coord),
                                                                                            ids = filter(pu.species.probabilities,species==i.species) %>% pull(pu)), # only occupied PUs
                                                  demand.points = dp,
                                                  species=as.integer(i.species))
    
}

# Attribute spaces for both species
attribute.spaces <- list()
attribute.spaces[[1]] <-
    AttributeSpaces(
        spaces = genetic_spaces,
        name = "genetic"
    )
# Set status of existing MR to 2 to lock them in (Marxan convention)
pus$status[pus$status==1] <- pus$status[pus$status==1]+1
pus$status <- as.integer(pus$status)

# Convert to PolySet and calculate boundary
pus %>%  convert2PolySet() -> polygons
boundary <- raptr::calcBoundaryData(polygons)

rap_data <- RapData(pu = st_drop_geometry(pus),
                    species = data.frame(name=c("Diplodus","Mullus")),
                    targets = targets,
                    pu.species.probabilities = pu.species.probabilities,
                    attribute.spaces = attribute.spaces,
                    boundary = boundary,
                    polygons = polygons)

# Define problem 
ro <- RapUnreliableOpts(BLM=0)
prob_gs <- RapUnsolved(ro, rap_data)
save(prob_gs, file="Results/Problem_definition_raptr.RData")
#maximum.targets(prob_gs)


################################################################################
# Evaluation of the current MR set for both species
################################################################################
load("Results/Problem_definition_raptr.RData")
selections <- which(prob_gs@data@pu$status==2)
prob_gs_current_MR <- update(prob_gs, b = selections)
# Space held
space.held(prob_gs_current_MR, y=NULL)
# Amount held
amount.held(prob_gs_current_MR, y=NULL)
# Total conservation cost
pus %>%
    st_drop_geometry %>%
    filter(status == 2) %>%
    pull(cost) %>%
    sum

# Space plot
## Add a column "selected" to pus to include the existing marine reserves ("selections")
pus %>% mutate(selected = 0) -> pus_with_selections
pus_with_selections$selected[selections] <- 1
png("Figures/space-plots.png",width=20,height=10,units="cm",res=300)
par(mfrow=c(1,2))
# Diplodus
## Coordinates of demand points of the first two PCA axes
prob_gs_current_MR@data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[,c(1,2)] %>%
    plot(pch=16, col="gray", main="Diplodus sargus",cex=0.5)
## Subset the planning units that are (i) existing marine reserves AND (ii) occupied by the species
selections_species <- which(filter(pus_with_selections,Diplodus_sargus == 1)$selected == 1)
## Plot them in green
prob_gs_current_MR@data@attribute.spaces[[1]]@spaces[[1]]@planning.unit.points@coords[selections_species,c(1,2)] %>%
    points(pch=16, col="green",cex=0.5)
# Mullus
## Coordinates of demand points of the first two PCA axes
prob_gs_current_MR@data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[,c(1,2)] %>%
    plot(pch=16, col="gray", main="Mullus surmuletus",cex=0.5)
## Subset the planning units that are (i) existing marine reserves AND (ii) occupied by the species
selections_species <- which(filter(pus_with_selections,Mullus_surmuletus == 1)$selected == 1)
## Plot them in green
prob_gs_current_MR@data@attribute.spaces[[1]]@spaces[[2]]@planning.unit.points@coords[selections_species,c(1,2)] %>%
    points(pch=16, col="green",cex=0.5)
dev.off()


################################################################################
# Extension of the current MR set: amount targets (with prioritizr)
################################################################################
prob_so <- problem(pus,
                   features = c("Diplodus_sargus","Mullus_surmuletus"),
                   cost_column = "cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(0.15) %>%
    add_binary_decisions() %>%
    add_locked_in_constraints(which(pus$status==2)) %>%
    add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
    add_gurobi_solver(gap=0, threads=threads)
res_so <- solve(prob_so, run_checks=T)

# Space held
space_held_Diplodus <- space_held_Mullus <- vector()
for (i.sol in 1 : 100) {
    if (i.sol %% 5 == 0) {
        cat(i.sol,"\n");
        flush.console()
    }
    selections <- which(res_so %>% pull(paste0("solution_",i.sol))==1)
    res_updated <- update(prob_gs, b = selections)
    space_held_Diplodus[i.sol] <- space.held(res_updated, y=1, species=1) %>% as.vector()
    space_held_Mullus[i.sol] <- space.held(res_updated, y=1, species=2) %>% as.vector()
}
space_held_Diplodus %>% summary
space_held_Mullus %>% summary

# Conservation cost
cost_solution <- vector()
for (i.sol in 1 : 100) {
    if (i.sol %% 25 == 0) {
        cat(i.sol,"\n");
        flush.console()
    }
    prob_so %>%
        eval_cost_summary(select(res_so,paste0("solution_",i.sol))) %>% pull(cost) ->
        cost_solution[i.sol]
}
summary(cost_solution)

save(prob_so,
     res_so,
     space_held_Diplodus,
     space_held_Mullus,
     cost_solution,
     file="Results/Results_Extension_Amount.RData")


