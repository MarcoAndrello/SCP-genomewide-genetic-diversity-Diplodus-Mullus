# Extend the current MR network
# Using quantile, 6 classes

rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(rgeoda)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tmap)

library(prioritizr)

# Number of assigned threads
threads <- 1L #as.integer(max(1, detectCores() - 2))

# Source functions
source(paste0(getwd(),"/functions/split.taxon.R"))

# Define classification method and number of classes
class_method <- "quantile"
num_classes <- 6

# Load Planning units
load("Planning_units.RData")

# Rescale mean cost to avoid numerical issues
pus$mean_cost <- pus$mean_cost / 1000

# Diplodus sargus
g_rast <- rast(paste0(getwd(),"/Results May_2023/Diplodus_allAxes_8068.grd"))
g_rast_values <- terra::extract(g_rast,pus_centroid)
pus %>% 
    select(Diplodus_sargus, Mullus_surmuletus, cost, status_0, mean_cost) -> 
    pus_split.taxon
for (i.axis in 1 : 17) {
    st_matrix <- split.taxon(x = g_rast_values[,i.axis],
                             num_classes = num_classes,
                             class_method = class_method,
                             rij_taxon = pull(pus,"Diplodus_sargus"),
                             name_feat = paste0("Diplodus_a",sprintf("%02d",i.axis)))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
}
# Mullus surmuletus
g_rast <- rast(paste0(getwd(),"/Results May_2023/Mullus_allAxes_2753.grd"))
g_rast_values <- terra::extract(g_rast,pus_centroid)
for (i.axis in 1 : 26) {
    st_matrix <- split.taxon(x = g_rast_values[,i.axis],
                             num_classes = num_classes,
                             class_method = class_method,
                             rij_taxon = pull(pus,"Mullus_surmuletus"),
                             name_feat = paste0("Mullus_a",sprintf("%02d",i.axis)))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
}
# Feature names
features <- names(pus_split.taxon)[c(grep("Diplodus_a",names(pus_split.taxon)),
                                     grep("Mullus_a",names(pus_split.taxon)))]

# Define and solve problem
prob <- problem(pus_split.taxon,
                features = features,
                cost_column = "mean_cost") %>%
                # cost_column = "cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(0.15) %>%
    add_binary_decisions() %>%
    add_locked_in_constraints(which(pus_split.taxon$status_0==1)) %>%
    # add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
    add_gurobi_solver(gap=0, threads=threads)
res <- solve(prob, run_checks=T)

# Evaluate targets, n and cost
eval_target_coverage_summary(prob,res["solution_1"]) %>%
    filter(relative_target > 0) %>%
    pull(relative_held) %>% hist
eval_n_summary(prob,res["solution_1"])
eval_cost_summary(prob,res["solution_1"])

# Plot solution
# Convert to factor for better plotting
res$solution_1_fac <- res$solution_1
res$solution_1_fac[which(is.na(res$solution_1_fac))] <- 0
res$solution_1_fac[which(res$solution_1_fac==0)] <- "not selected"
res$solution_1_fac[which(res$solution_1_fac==1)] <- "selected"
res$solution_1_fac[which(res$status_0==1)] <- "existing MR"
res$solution_1_fac <- factor(res$solution_1_fac,levels=,c("existing MR","not selected","selected"))
# res$status_0_fac <- res$status_0
# res$status_0_fac[which(res$status_0_fac == 0)] <- NA
# res$status_0_fac <- factor(res$status_0_fac)
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(res)) -> countries
png("Map_prioritizr_twospecies.png",width=20,height=12,res=300,units="cm")
tm_shape(res) +
    tm_fill("solution_1_fac", palette=c("blue","gray","red"), legend.is.portrait = F) +
    tm_legend(legend.outside = T, legend.outside.position = "bottom", legend.titl) +
    tm_shape(countries, bbox = res) +
    tm_polygons(col="lightgray")
dev.off()

### why locked_in PUs are not in the solution ???
# res %>% filter(status_0 == 1) %>% select(Mullus_surmuletus, mean_cost, status_0, solution_1) %>% print(n=500)
# res %>% filter(status_0 == 1) %>% select(Diplodus_sargus, mean_cost, status_0, solution_1) %>% boxplot(mean_cost ~ solution_1,data=.)
# res %>% st_drop_geometry %>% select(paste0("solution_",c(1:100))) %>% rowSums(na.rm=T) -> res$selection_frequency
# ne_countries(returnclass = "sf") %>% st_transform(st_crs(res)) -> countries
# ne_coastline(scale = 50, returnclass = "sf") %>% st_transform(st_crs(res)) -> coastline
# plot(st_geometry(coastline),extent=res)
# plot(res["selection_frequency"], border=NA, add=T)

# # Replacement importance
# replacement_importance <- eval_replacement_importance(prob,res["solution_1"], run_checks=F)
# save(replacement_importance,file=paste0("Replacement_importance_twospecies.RData"))


# Comparison with speices-only SCP
prob_so <- problem(pus,
                features = c("Diplodus_sargus","Mullus_surmuletus"),
                cost_column = "mean_cost") %>%
    add_min_set_objective() %>%
    add_relative_targets(0.15) %>%
    add_binary_decisions() %>%
    add_locked_in_constraints(which(pus_split.taxon$status_0==1)) %>%
    add_gurobi_solver(gap=0, threads=threads)
res_so <- solve(prob_so, run_checks=T)

eval_target_coverage_summary(prob_so,res_so["solution_1"]) %>%
    pull(relative_held)
eval_n_summary(prob_so,res_so["solution_1"])
eval_cost_summary(prob_so,res_so["solution_1"])

# Does species-only SCP reach targets for genetic SCP?
eval_target_coverage_summary(prob,res_so["solution_1"]) %>%
    pull(relative_held) %>%
    sort
