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
library(tmap)
library(RColorBrewer)

library(prioritizr)

# Number of assigned threads
threads <- 1L #as.integer(max(1, detectCores() - 2))

# Source functions
source(paste0(getwd(),"/functions/split.taxon.R"))

# Load Planning units
load("Planning_units.RData")

# Rescale mean cost to avoid numerical issues
pus$mean_cost <- pus$mean_cost / 1000

g_rast <- rast(paste0(getwd(),"/Results May_2023/Diplodus_allAxes_8068.grd"))
g_rast_values_Diplodus <- terra::extract(g_rast,pus_centroid)
g_rast <- rast(paste0(getwd(),"/Results May_2023/Mullus_allAxes_2753.grd"))
g_rast_values_Mullus <- terra::extract(g_rast,pus_centroid)
rm(g_rast)

# Single PCA axis
prob <- res <- list()
v.class_method <- c("quantile","equal","natural","stddev")
v.num_classes <- c(3,6,12)
i.prob <- 1
# Loop on classification method and number of classes
for (i.class_method in 1 : length(v.class_method)) {
    class_method <- v.class_method[i.class_method]
    for (i.num_classes in  1 : length(v.num_classes)) {
        num_classes <- v.num_classes[i.num_classes]
        if (class_method == "stddev" & num_classes != 6) next # only 6 classes for stddev
        # Diplodus sargus
        pus %>% 
            select(Diplodus_sargus, Mullus_surmuletus, cost, status_0, mean_cost) -> 
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
        # Feature names
        features <- names(pus_split.taxon)[c(grep("Diplodus_a",names(pus_split.taxon)),
                                             grep("Mullus_a",names(pus_split.taxon)))]
        
        # Define and solve problem
        prob[[i.prob]] <- problem(pus_split.taxon,
                                  features = features,
                                  cost_column = "mean_cost") %>%
            # cost_column = "cost") %>%
            add_min_set_objective() %>%
            add_relative_targets(0.15) %>%
            add_binary_decisions() %>%
            add_locked_in_constraints(which(pus_split.taxon$status_0==1)) %>%
            # add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
            add_gurobi_solver(gap=0, threads=threads)
        res[[i.prob]] <- solve(prob[[i.prob]], run_checks=T)
        i.prob <- i.prob + 1
    }
}

# Multiple PCA axes
prob_multi <- res_multi <- list()
i.num_classes <- 1
for (i.num_classes in 1 : 2){
    pus %>% 
        select(Diplodus_sargus, Mullus_surmuletus, cost, status_0, mean_cost) -> 
        pus_split.taxon
    # Diplodus
    num_classes <- c(2,5)[i.num_classes]
    st_matrix <- split.taxon.multi(x = g_rast_values_Diplodus[,1:17],
                                   num_classes = num_classes,
                                   rij_taxon = pull(pus,"Diplodus_sargus"),
                                   name_feat = paste0("Diplodus_m17"))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
    # Mullus
    num_classes <- c(3,7)[i.num_classes]
    st_matrix <- split.taxon.multi(x = g_rast_values_Mullus[,1:26],
                                   num_classes = num_classes,
                                   rij_taxon = pull(pus,"Mullus_surmuletus"),
                                   name_feat = paste0("Mullus_m26"))
    pus_split.taxon <- cbind(pus_split.taxon, st_matrix)
    # Feature names
    features <- names(pus_split.taxon)[c(grep("Diplodus_m",names(pus_split.taxon)),
                                         grep("Mullus_m",names(pus_split.taxon)))]
    ## Solve problem
    prob_multi[[i.num_classes]] <- problem(pus_split.taxon,
                 features = features,
                 cost_column = "mean_cost") %>%
        add_min_set_objective() %>%
        add_relative_targets(0.15) %>%
        add_binary_decisions() %>%
        # add_gap_portfolio(number_solutions = 100, pool_gap = 0.02) %>%
        add_gurobi_solver(gap=0, threads=threads)
    res_multi[[i.num_classes]] <- solve(prob_multi[[i.num_classes]], run_checks=F)
}

prob <- c(prob, prob_multi)
res <- c(res, res_multi)
names(prob) <- names(res) <- 
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

# Evaluate targets, n and cost
n <- cost <- vector()
for (i.prob in 1 : length(prob)) {
    n[i.prob] <- 
        eval_n_summary(prob[[i.prob]],res[[i.prob]]["solution_1"]) %>% pull(n)
    cost[i.prob] <- 
        eval_cost_summary(prob[[i.prob]],res[[i.prob]]["solution_1"]) %>% pull(cost)
}

res[[1]] %>% select(Diplodus_sargus, Mullus_surmuletus, status_0, mean_cost, solution_1) -> results
res[[2]]$solution_1 -> results$solution_2
res[[3]]$solution_1 -> results$solution_3
res[[4]]$solution_1 -> results$solution_4
res[[5]]$solution_1 -> results$solution_5
res[[6]]$solution_1 -> results$solution_6
res[[7]]$solution_1 -> results$solution_7
res[[8]]$solution_1 -> results$solution_8
res[[9]]$solution_1 -> results$solution_9
res[[10]]$solution_1 -> results$solution_10
res[[11]]$solution_1 -> results$solution_11
res[[12]]$solution_1 -> results$solution_12
results %>% mutate(sel_frequency =
                       (solution_1 +
                       solution_2 +
                       solution_3 +
                       solution_4 +
                       solution_5 +
                       solution_6 +
                       solution_7 +
                       solution_8 +
                       solution_9 +
                       solution_10 +
                       solution_11 +
                       solution_12) / 12) -> results

# Assign dummy value to locked_in PUs
results$sel_frequency %>% cut(breaks=seq(0,1,0.2),include.lowest=T,right=F) -> results$sel_frequency_fac                      
summary(results$sel_frequency_fac)
results$selection_frequency <- factor(results$sel_frequency_fac,
                                      levels = c("Existing",levels(results$sel_frequency_fac)))
results$selection_frequency[results$status_0 == 1] <- "Existing"
summary(results$selection_frequency)

ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(results)) -> countries
png("Map_prioritizr_twospecies.png",width=20,height=12,res=300,units="cm")
tm_shape(results) +
    tm_fill("selection_frequency", palette=c("gold",brewer.pal(5,"Greens")), legend.is.portrait = F) +
    tm_legend(legend.outside = T, legend.outside.position = "bottom") +
    tm_shape(countries, bbox = res) +
    tm_polygons(col="lightgray")
dev.off()

n / 551
cost / 1047296
par(mar=c(5,4,2,2))
barplot(n, names.arg=names(prob),las=2,cex.names=0.75,cex.axis=0.75,
        main="Number of PUs", ylab="Number")
abline(h=551,lty=2)
barplot(cost, names.arg=names(prob),las=2,cex.names=0.75,cex.axis=0.75,
        main="Conservation cost", ylab="Thousand of euros")
abline(h=1047296,lty=2)

# Comparison with species-only SCP
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
n_met <- percent_met <- vector()
for (i.prob in 1 : 12) {
    eval_target_coverage_summary(prob[[i.prob]],res_so["solution_1"]) %>%
    pull(met) %>% which %>% length -> n_met[i.prob]
    eval_target_coverage_summary(prob[[i.prob]],res_so["solution_1"]) %>%
        filter(relative_target > 0) %>% nrow -> n_cf
    percent_met[i.prob] <- n_met[i.prob] / n_cf
}

    
