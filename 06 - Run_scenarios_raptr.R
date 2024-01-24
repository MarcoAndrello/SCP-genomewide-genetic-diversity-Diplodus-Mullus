# Run scenarios prioritizr
rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(raptr)
library(parallel)
library(tictoc)

load("Results/Problem_definition_raptr.RData")


################################################################################
# Extension of the current MR set: amount + space targets (with raptr)
################################################################################


tic()
res_gs <- solve(prob_gs, Threads = 1L, verbose=T, NumericFocus= 3L,
                MIPGap=0.02, NumberSolutions=1L)
toc()
# runs for 7 hours then throws an out-of-memory error
#space.held(res_gs, y=NULL)
#amount.held(res_gs, y=NULL)




################################################################################
# Approximating the genetic spaces with a lower number of demand points
################################################################################

# Using 50% of the demand points of the gold standard

set.seed(20231214)

# Create a new rapdata object from that of the golden standard
rap_data_50 <- prob_gs@data
prob_50gs <- res_50gs <- list()
# Loop on replicates: warning: takes 2-3 hours per replicate
for (i.perm in 1 : 20) {
    # Random sample of demand points for each species
    id_Diplodus <- sample(2253,2253/2)
    id_Mullus <- sample(3613,3613/2)
    # Populate the new RapData object with the coordinates and weights of the new demand points
    ## For Diplodus
    rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
    rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <-
        rep(1/length(id_Diplodus),length(id_Diplodus))
    ## For Mullus
    rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
    rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <-
        rep(1/length(id_Mullus),length(id_Mullus))
    # Define raptr problem and solve it
    ro <- RapUnreliableOpts(BLM=0)
    prob_50gs[[i.perm]] <- RapUnsolved(ro, rap_data_50)
    tic()
    res_50gs[[i.perm]] <- solve(prob_50gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                                MIPGap=0.02, NumberSolutions=1L)
    toc()
}
save(prob_50gs, file=paste0("Results/Problems_raptr_50gs.RData"))
save(res_50gs, file=paste0("Results/Results_raptr_50gs.RData"))


# Using 20% of the demand points of the gold standard

set.seed(20231214)

# Create a new rapdata object from that of the golden standard
rap_data_20 <- prob_gs@data
prob_20gs <- res_20gs <- list()
# Loop on replicates: warning: takes 2-3 hours per replicate
for (i.perm in 1 : 20) {
    # Random sample of demand points for each species
    id_Diplodus <- sample(2253,2253*0.2)
    id_Mullus <- sample(3613,3613*0.2)
    # Populate the new RapData object with the coordinates and weights of the new demand points
    ## For Diplodus
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <-
        rep(1/length(id_Diplodus),length(id_Diplodus))
    ## For Mullus
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <-
        rep(1/length(id_Mullus),length(id_Mullus))
    # Define raptr problem and solve it
    ro <- RapUnreliableOpts(BLM=0)
    prob_20gs[[i.perm]] <- RapUnsolved(ro, rap_data_20)
    tic()
    res_20gs[[i.perm]] <- solve(prob_20gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                                MIPGap=0.02, NumberSolutions=1L)
    toc()
}
save(prob_20gs, file=paste0("Results/Problems_raptr_20gs.RData"))
save(res_20gs, file=paste0("Results/Results_raptr_20gs.RData"))
