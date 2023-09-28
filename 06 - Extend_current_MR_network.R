# Extend the current MR network
# Using equal interval, 12 classes

rm(list=ls())

species <- "Diplodus"

library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(rgeoda)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

#############################################################################
# Prioritizr
#############################################################################
library(prioritizr)

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

class_method <- "quantile"
num_classes <- 6

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

# Feature names
features <- names(pus_split.taxon)[grep(name_species_feat,names(pus_split.taxon))]
# rescale cost to avoid numerical issues
pus_split.taxon$mean_cost <- pus_split.taxon$mean_cost / 1000

# Define and solve problem
prob <- problem(pus_split.taxon,
                features = features,
                cost_column = "mean_cost") %>%
                # cost_column = "cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(0.15) %>%
    add_binary_decisions() %>%
    add_locked_in_constraints(which(pus_split.taxon$status_0==1)) %>%
    add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
    add_gurobi_solver(gap=0.02, threads=threads)
res <- solve(prob, run_checks=F)
plot(res["solution_1"],border=NA)

res %>% filter(status_0 == 1) %>% select(Diplodus_sargus, mean_cost, status_0, solution_1) %>% print(n=500)
res %>% filter(status_0 == 1) %>% select(Diplodus_sargus, mean_cost, status_0, solution_1) %>% boxplot(mean_cost ~ solution_1,data=.)

### why locked_in PUs are not in the solution ???
res %>% st_drop_geometry %>% select(paste0("solution_",c(1:100))) %>% rowSums(na.rm=T) -> res$selection_frequency
ne_countries(returnclass = "sf") %>% st_transform(st_crs(res)) -> countries
ne_coastline(scale = 50, returnclass = "sf") %>% st_transform(st_crs(res)) -> coastline
plot(st_geometry(coastline),extent=res)
plot(res["selection_frequency"], border=NA, add=T)

replacement_importance <- eval_replacement_importance(prob,res["solution_1"], run_checks=F)
save(replacement_importance,file=paste0("Replacement_importance_prioritizr_",species,".RData"))
plot(st_geometry(coastline),extent=res)
plot(replacement_importance["rc"], border=NA, add=T)





#####################################################################################
# Prioritizr, targets met 
#####################################################################################
load(paste0("Results_prioritizr_",species,".RData"))
rm(results)

# For each scenario, I analyze whether the other scenarios meet its target,by calculating relative targets over conservation features
amount_held <- array(NA,length(problems))
i.prob <- 1
for (i.prob in 1 : length(problems)) {
    eval_target_coverage_summary(problems[[i.prob]],
                                 select(pus, status_0)) %>%
        filter(relative_target > 0) %>%
        pull(relative_held) %>%
        mean() ->
        amount_held[i.prob]
}
amount_held <- tibble(problem=names(problems),amount_held=amount_held)
amount_held$problem <- factor(amount_held$problem,
                              levels=names(problems))

png(paste0("Evaluate_MR_prioritizr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)
theme_set(theme_classic())
ggplot(amount_held,aes(x=problem,y=amount_held)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    geom_hline(yintercept=0.15, linetype="dotted") +
    ylim(0,0.15)
dev.off()
rm(problems, amount_held, i.prob)
#    
    
    
#####################################################################################
# Raptr
#####################################################################################
load(paste0("Results_raptr_",species,".RData"))

selections <- which(pus$status_0 == 1)
space_held <- distance_to_maximum <- vector()
for (i.prob in 1 : length(results)){
    cat(i.prob,"\n"); flush.console()
    maximum_targets_problem <- maximum.targets(results[[i.prob]])$proportion
    res_updated <- update(results[[i.prob]], b = selections)
    space_held_feature <- space.held(res_updated, y=1, space=NULL) %>% as.vector()
    space_held[i.prob] <- mean(space_held_feature)
    # distance_to_maximum[i.prob] <- mean(maximum_targets_problem - space_held_feature)
}
rm(maximum_targets_problem, space_held_feature, selections, res_updated)
space_held <- tibble(problem=names(problems),space_held=space_held)
space_held$problem <- factor(space_held$problem,
                            levels=names(problems))
png(paste0("Evaluate_MR_raptr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)
theme_set(theme_classic())
ggplot(space_held ,aes(x=problem,y=space_held )) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) 
dev.off()

# distance_to_maximum <- tibble(problem=names(problems),distance_to_maximum=distance_to_maximum)
# distance_to_maximum$problem <- factor(distance_to_maximum$problem,
#                               levels=names(problems))
# png(paste0("Evaluate_MR_raptr_distance_",species,".png"),width=7.5,height=7.5,units="cm",res=300)
# theme_set(theme_classic())
# ggplot(distance_to_maximum ,aes(x=problem,y=distance_to_maximum )) +
#     geom_point() +
#     theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) 
# dev.off()
rm(problems, results, space_held, distance_to_maximum, i.prob)

