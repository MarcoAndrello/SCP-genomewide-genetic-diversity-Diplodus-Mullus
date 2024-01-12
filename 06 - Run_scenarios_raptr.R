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
rap_data_50 <- prob_gs@data
prob_50gs <- res_50gs <- list()
for (i.perm in 16 : 20) {
    id_Diplodus <- sample(2253,2253/2)
    id_Mullus <- sample(3613,3613/2)
    if(i.perm < 6) next
    rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
    rap_data_50@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <-
        rep(1/length(id_Diplodus),length(id_Diplodus))
    rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
    rap_data_50@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <-
        rep(1/length(id_Mullus),length(id_Mullus))
    
    ro <- RapUnreliableOpts(BLM=0)
    prob_50gs[[i.perm]] <- RapUnsolved(ro, rap_data_50)
    tic()
    res_50gs[[i.perm]] <- solve(prob_50gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                                MIPGap=0.02, NumberSolutions=1L)
    toc()
    save(prob_50gs, res_50gs, file=paste0("Results/Results_raptr_50gs_perm",i.perm,".RData"))
}
# save(prob_50gs, res_50gs, file=paste0("Results/Results_raptr_50gs.RData"))

# Using 20% of the demand points of the gold standard
set.seed(20231214)
rap_data_20 <- prob_gs@data
prob_20gs <- res_20gs <- list()
for (i.perm in 1 : 20) {
    # Sample 20% of demand points randomly
    id_Diplodus <- sample(2253,2253*0.2)
    id_Mullus <- sample(3613,3613*0.2)
    if(i.perm < 13) next
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[1]]@demand.points@coords[id_Diplodus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[1]]@demand.points@weights <-
        rep(1/length(id_Diplodus),length(id_Diplodus))
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords <-
        prob_gs@data@attribute.spaces[[1]]@spaces[[2]]@demand.points@coords[id_Mullus,]
    rap_data_20@attribute.spaces[[1]]@spaces[[2]]@demand.points@weights <-
        rep(1/length(id_Mullus),length(id_Mullus))
    # Define problem and solve it
    ro <- RapUnreliableOpts(BLM=0)
    prob_20gs[[i.perm]] <- RapUnsolved(ro, rap_data_20)
    tic()
    res_20gs[[i.perm]] <- solve(prob_20gs[[i.perm]], Threads = 1L, verbose=T, NumericFocus= 3L,
                                MIPGap=0.02, NumberSolutions=1L)
    toc()
    save(prob_20gs, res_20gs, file=paste0("Results/Results_raptr_20gs_perm",i.perm,".RData"))
}
# save(prob_20gs, res_20gs, file="Results/Results_raptr_20gs.RData")

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

load("Results/Results_raptr_50gs.RData")
load("Results/Results_raptr_20gs.RData")

problems_raptr <- c(prob_50gs,prob_20gs)
results_raptr <- c(res_50gs,res_20gs)
names(problems_raptr) <- names(results_raptr) <- 
    c("raptr_50perc_repl1",
      "raptr_50perc_repl2",
      "raptr_50perc_repl3",
      "raptr_50perc_repl4",
      "raptr_50perc_repl5",
      "raptr_20perc_repl1",
      "raptr_20perc_repl2",
      "raptr_20perc_repl3",
      "raptr_20perc_repl4",
      "raptr_20perc_repl5")

# Save
save(problems_raptr,
     results_raptr,
     file=paste0(getwd(),"/Results/Results_raptr.RData"))

