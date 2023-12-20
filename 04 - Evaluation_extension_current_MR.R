# Evaluation and extension of the current set of MR
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(raptr)
library(prioritizr)
library(parallel)
library(tictoc)

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
    names(pus_g_values)[1:ncol(g_rast_values)] <- paste0("spca",sprintf("%02d",1:ncol(g_rast_values)))
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
maximum.targets(prob_gs)


################################################################################
# Evaluation of the current MR set for both species
################################################################################

selections <- which(rap_data@pu$status==2)
prob_gs_current_MR <- update(prob_gs, b = selections)
space.held(prob_gs_current_MR , y=NULL)
amount.held(prob_gs_current_MR , y=NULL)
pus %>% st_drop_geometry %>% filter(status == 2) %>% pull(cost) %>% sum



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



################################################################################
# Extension of the current MR set: amount + space targets (with raptr)
################################################################################
tic()
res_gs <- solve(prob_gs, Threads = 1L, verbose=T, NumericFocus= 3L,
                MIPGap=0.02, NumberSolutions=1L)
toc()
space.held(res_gs, y=NULL)
amount.held(res_gs, y=NULL)
# spp.plot(res_gs, species=1, y=1) #, pu.color.palette = c("grey30", "green", "red", "red"))

save(prob_gs,
     res_gs,
     file=paste0(getwd(),"/Results/Results_raptr_100PORTFOLIO.RData"))


################################################################################
# Approximating the genetic spaces with a lower number of demand points
################################################################################
# Using 50% of the demand points of the gold standard
set.seed(20231214)
rap_data_50 <- rap_data
prob_50gs <- res_50gs <- list()
i.perm <- 2
id_Diplodus <- sample(2253,2253/2)
id_Mullus <- sample(3613,3613/2)
rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
    rap_data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <- rep(1/length(id_Diplodus),length(id_Diplodus))
rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
    rap_data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <- rep(1/length(id_Mullus),length(id_Mullus))

ro <- RapUnreliableOpts(BLM=0)
prob_50gs[[i.perm]] <- RapUnsolved(ro, rap_data_50)
tic()
res_50gs[[i.perm]] <- solve(prob_50gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                MIPGap=0.02, NumberSolutions=1L)
toc()
save(prob_50gs, res_50gs, file="Results/Results_raptr_50gs_perm2.RData")

# Using 20% of the demand points of the gold standard
set.seed(20231214)
rap_data_20 <- rap_data
prob_20gs <- res_20gs <- list()
for (i.perm in 1 : 5) {
    # Sample 20% of demand points randomly
    id_Diplodus <- sample(2253,2253*0.2)
    id_Mullus <- sample(3613,3613*0.2)
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
        rap_data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <- rep(1/length(id_Diplodus),length(id_Diplodus))
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
        rap_data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <- rep(1/length(id_Mullus),length(id_Mullus))
    # Define problem and solve it
    ro <- RapUnreliableOpts(BLM=0)
    prob_20gs[[i.perm]] <- RapUnsolved(ro, rap_data_20)
    tic()
    res_20gs[[i.perm]] <- solve(prob_20gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                                MIPGap=0.02, NumberSolutions=1L)
    toc()
}
save(prob_20gs, res_20gs, file="Results/Results_raptr_20gs.RData")

load("Results/Results_raptr_20gs.RData")
space_held_Diplodus <- space_held_Mullus <- vector()
for (i.perm in 1 : 5) {
    cat(i.perm,"\n")
    selections <- which(as.vector(res_20gs[[i.perm]]@results@selections) == 1)
    res_updated <- update(prob_gs, b = selections)
    space_held_Diplodus[i.perm] <- space.held(res_updated, y=1, species=1)
    space_held_Diplodus[i.perm] <- space.held(res_updated, y=1, species=2)
}
data.frame(space_held_Diplodus = space_held_Diplodus,
           space_held_Mullus = space_held_Mullus,
           solution = "raptr_20percent",
           sol = as.integer(1:5))



