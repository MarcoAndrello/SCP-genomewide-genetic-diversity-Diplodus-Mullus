# Assess scenarios raptr
# (1) pairwise distance between PU selection frequencies
# (2) Targets met by solutions found with different scenarios
# (3) Conservation cost
# (4) Maximum targets

rm(list=ls())

species <- "Diplodus"

library(tidyverse)
library(sf)
library(ade4)
library(vegan)
library(raptr)

load(paste0("Results_raptr_",species,".RData"))


#####################################################################################
# (1) pairwise distance between PU selection frequencies
#####################################################################################
# i) I calculate a selection frequency for each of the two scenarios
# ii) I use the 2 selection frequencies to calculate a Jaccard distance
# iii) I assess significance by permuting solution vectors between the 2 scenarios 
jaccard_distance <- pvalue <- array(NA,dim=c(length(results),length(results)))
dimnames(jaccard_distance) <- dimnames(pvalue) <- list("1" = names(results), "2" = names(results))

par(mar=c(1,1,4,1),mfrow=c(4,4))
i.res <- 1; j.res <- 2
for (i.res in 1 : (length(results)-1)) {
    for (j.res in (i.res+1) : length(results)) {
        cat(i.res,j.res,"\n"); flush.console()
        twogroups <- cbind(
            results[[i.res]]@results@selections %>% t(),
            results[[j.res]]@results@selections %>% t()
        )
        selection_frequency_1 <- rowSums(twogroups[,1:100]) / 100
        selection_frequency_2 <- rowSums(twogroups[,101:200]) / 100
        vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_observed
        # # Permute columns
        # num_perm <- 1000
        # distance_permuted <- rep(NA,num_perm)
        # for (i.perm in 1 : num_perm) {
        #     if (i.perm %% 100 == 0) {
        #         cat(i.perm,"of",num_perm,"\n")
        #         flush.console()
        #     }
        #     twogroups_permuted <- twogroups[,sample(200)]
        #     selection_frequency_1 <- rowSums(twogroups_permuted[,1:100]) / 100
        #     selection_frequency_2 <- rowSums(twogroups_permuted[,101:200]) / 100
        #     vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_permuted[i.perm]
        # }
        # hist(distance_permuted,
        #      main = names(results)[c(i.res,j.res)],
        #      xlim=c(min(distance_permuted),
        #             max((distance_observed+1), max(distance_permuted))))
        # abline(v=distance_observed,col="red")
        # jaccard_distance[i.res,j.res] <- distance_observed
        # pvalue[i.res,j.res] <- 1 - (ecdf(distance_permuted)(distance_observed))
    }
}
#save(jaccard_distance, pvalue, file=paste0("Jaccard_distance_raptr_",species,".RData"))
# Plot a tree showing distances between solutions
load(paste0("Jaccard_distance_raptr_",species,".RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png(paste0("hclust_raptr_",species,".png"),width=9,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main=paste("raptr",species),xlab="",sub="")
dev.off()

# For Diplodus
# Distance among single PCA axes solutions
jaccard_distance[1:10,1:10] %>% as.vector() %>% summary()
pvalue[1:10,1:10] %>% as.vector() %>% summary()
# Distance among combined PCA axes solutions
jaccard_distance[11:13,11:13] %>% as.vector() %>% summary()
pvalue[11:13,11:13] %>% as.vector() %>% summary()
# Distance between single PCA axes and combined PCA axes
jaccard_distance[1:10,11:13] %>% as.vector() %>% summary()
pvalue[1:10,11:13] %>% as.vector() %>% summary()

# For Mullus
# Distance among solutions except single_natural_12
jaccard_distance[-9,-9] %>% as.vector() %>% summary()
pvalue[-9,-9] %>% as.vector() %>% summary()
load("Planning_units.RData")
## (rerun selection_frequency_1 etc for selected scenarios)
pus$diff_sel <- selection_frequency_1 - selection_frequency_2
plot(pus["diff_sel"],border=NA,
     breaks=seq(-1,1,0.25),
     pal=c("darkred","red","orange","gray","gray","darkgreen","green","lightgreen"))

#####################################################################################
# (2) Space held by solutions found with different scenarios
#####################################################################################


# For each solution, I calculate the space held in other planning problems
list_space_held <- list()
plots <- list()
theme_set(theme_classic())
i.prob <- j.prob <- 1
for (i.prob in 1 : length(results)) {
    k <- 1
    for (j.prob in 1 : length(problems)) {
        for (i.sol in 1 : 100) {
            if (i.sol %% 25 == 0) {
                cat(i.prob,j.prob,i.sol,"\n");
                flush.console()
            }
            if(i.prob == j.prob) {
                space_held_solution <- space.held(results[[i.prob]], y=i.sol,space=NULL) %>% as.vector()
            } else {
                selections <- which(results[[i.prob]]@results@selections[i.sol,]==1)
                res_updated <- update(results[[j.prob]], b = selections)
                space_held_solution <- space.held(res_updated, y=1,space=NULL) %>% as.vector()
            }
            data.frame(space_held = space_held_solution,
                       problem = names(problems)[[j.prob]],
                       solution = names(results)[[i.prob]],
                       sol = as.integer(i.sol)) ->
                list_space_held[[k]]
            k <- k + 1
        }
    }
    space_held <- bind_rows(list_space_held)
    space_held$problem <- factor(space_held$problem, levels=names(problems))
    space_held$solution <- factor(space_held$solution, levels=names(results))
    
    ggplot(space_held) +
        geom_violin(aes(x=problem, y=space_held), fill="black") +
        theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
        ggtitle(paste(species,"- solutions of",names(problems)[[i.prob]])) ->
        plots[[i.prob]]
}
save(plots,file=paste0("Plots_space_held_",species,".RData"))
load(paste0("Plots_space_held_",species,".RData"))
png(paste0("Space_held_raptr_",species,".png"),width=15,height=10,units="cm",res=300)
plots[[2]] + geom_hline(yintercept=0.75,linetype="dotted")
dev.off()

# # For each scenario, I calculate the space held under different solutions
# list_space_held <- list()
# plots <- list()
# theme_set(theme_classic())
# i.prob <- j.prob <- 1
# for (i.prob in 1 : length(problems)) {
#     k <- 1
#     for (j.prob in 1 : length(results)) {
#         for (i.sol in 1 : 100) {
#             if (i.sol %% 25 == 0) {
#                 cat(i.prob,j.prob,i.sol,"\n");
#                 flush.console()
#             }
#             if(i.prob == j.prob) {
#                 space_held_solution <- space.held(results[[i.prob]], y=i.sol,space=NULL) %>% as.vector()
#             } else {
#                 selections <- which(results[[j.prob]]@results@selections[i.sol,]==1)
#                 res_updated <- update(results[[i.prob]], b = selections)
#                 space_held_solution <- space.held(res_updated, y=1,space=NULL) %>% as.vector()
#             }
#             data.frame(space_held = space_held_solution,
#                        problem = names(problems)[[i.prob]],
#                        solution = names(results)[[j.prob]],
#                        sol = as.integer(i.sol)) ->
#                 list_space_held[[k]]
#             k <- k + 1
#         }
#     }
#     space_held <- bind_rows(list_space_held)
#     space_held$problem <- factor(space_held$problem, levels=names(problems))
#     space_held$solution <- factor(space_held$solution, levels=names(results))
#     ggplot(space_held) +
#         geom_violin(aes(x=solution, y=space_held), fill="black") +
#         theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
#         ggtitle(paste(species,names(problems)[[i.prob]])) ->
#         plots[[i.prob]]
# }
# save(plots,file="plots_raptr_mullus.RData")



# space_held <- distance_to_maximum <- array(NA,c(length(problems),length(results),100))
# i.prob <- j.prob <- 1
# for (i.prob in 1 : length(problems)) {
#     maximum_targets_problem <- maximum.targets(results[[i.prob]])$proportion
#     for (j.prob in 1 : length(results)) {
#         for (i.sol in 1 : 100) {
#             if (i.sol %% 25 == 0) {
#                 cat(i.prob,j.prob,i.sol,"\n");
#                 flush.console()
#             }
#             if(i.prob == j.prob) {
#                 space_held_solution <- space.held(results[[i.prob]], y=i.sol,space=NULL) %>% as.vector()
#             } else {
#                 selections <- which(results[[j.prob]]@results@selections[i.sol,]==1)
#                 res_updated <- update(results[[i.prob]], b = selections)
#                 space_held_solution <- space.held(res_updated, y=1,space=NULL) %>% as.vector()
#             }
#             space_held[i.prob,j.prob,i.sol] <- mean(space_held_solution)
#             distance_to_maximum[i.prob,j.prob,i.sol] <- mean(maximum_targets_problem - space_held_solution)
#         }
#     }
# }
# save(space_held,distance_to_maximum,file=paste0("Space_held_raptr_",species,".RData"))
# load(paste0("Space_held_raptr_",species,".RData"))
# # Reshaping space_held and distance_to_maximum into a dataframe
# space_held_list <- distance_to_maximum_list <- list()
# for(i.prob in 1 : 13) {
#     space_held_list[[i.prob]] <- 
#         space_held[i.prob,,] %>%
#         as.data.frame() %>%
#         mutate(solution=factor(names(results),levels=names(results))) %>% 
#         pivot_longer(cols=1:100) %>%
#         mutate(space_held = value, problem = names(problems)[i.prob]) %>%
#         select(-value)
#     # distance_to_maximum_list[[i.prob]] <-
#     #     distance_to_maximum[i.prob,,] %>%
#     #     as.data.frame() %>%
#     #     mutate(scenario=factor(names(results),levels=names(results))) %>%
#     #     pivot_longer(cols=1:100) %>%
#     #     mutate(distance_to_maximum = value, problem = names(problems)[i.prob]) %>%
#     #     select(-value)
# }
# space_held <- bind_rows(space_held_list)
# space_held$problem <- factor(space_held$problem,
#                             levels=names(problems))
# rm(space_held_list)
# # distance_to_maximum <- bind_rows(distance_to_maximum_list)
# # distance_to_maximum$problem <- factor(distance_to_maximum$problem,
# #                                       levels=names(problems))
# # rm(distance_to_maximum_list)
# 
# # Summarise by problem and scenario
# space_held %>% group_by(solution, problem) %>% summarise(mean_space_held=mean(space_held)) %>% 
#     mutate(space_held=mean_space_held) %>% select(-mean_space_held) -> mean_space_held
# ## Define method and num_clasees column
# mean_space_held$method = c(rep("quantile",13*3),
#                             rep("equal",13*3),
#                             rep("natural",13*3),
#                             rep("sd",13*1),
#                             rep("hypervolume",13*3))
# mean_space_held$method = factor(mean_space_held$method,
#                                  levels=c("quantile","equal","natural","sd","hypervolume"))
# mean_space_held$num_classes=c(rep(
#     rep(c("low","medium","high"),each=13),
#     3),
#     rep("medium",13),
#     rep("low",13),
#     rep("medium",13),
#     rep("high",13)
# )
# mean_space_held$num_classes = factor(mean_space_held$num_classes,levels=c("low","medium","high"))
# # Plot
# png(paste0("space_held_raptr_",species,".png"),width=20,height=11,units="cm",res=300)
# theme_set(theme_classic())
# ggplot(mean_space_held,aes(x=problem,y=space_held)) +
#     geom_point(aes(colour=method, shape=num_classes),position=position_dodge(width=0.5)) +
#     theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
#     scale_fill_brewer(type="qual") #+
# dev.off()


#####################################################################################
# (3) Conservation cost
#####################################################################################
cost_list <- list()
for (i.prob in 1 : length(results)) {
    cost_list[[i.prob]] <- results[[i.prob]]@results@summary$Cost
}
names(cost_list) <- names(results)
bind_cols(cost_list) %>%
    mutate(solution = c(1:100)) %>%
    pivot_longer(col=names(problems),
                 names_to="problem",
                 values_to="cost",
                 cols_vary="slowest") ->
    cost
cost$problem <- factor(cost$problem,
                       levels=names(problems))
theme_set(theme_classic())
png(paste0("Cost_raptr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)  
ggplot(cost,aes(x=problem,y=cost)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    ggtitle(paste0(paste("raptr",species)))
dev.off()


#####################################################################################
# (4) Maximum targets
#####################################################################################
maximum_targets_list <- lapply(problems, maximum.targets)
maximum_targets <- bind_rows(maximum_targets_list, .id="problem")
maximum_targets$problem <- factor(maximum_targets$problem,
                                  levels=names(problems))
save(maximum_targets,file=paste0("Maximum_targets_raptr_",species,".RData"))
load(paste0("Maximum_targets_raptr_",species,".RData"))
# Plot
png(paste0("Maximum_target_raptr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)    
theme_set(theme_classic())
ggplot(maximum_targets,aes(x=problem,y=proportion)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    ggtitle(paste0(paste("raptr",species)))
dev.off()

# For Diplodus: show why some maximum targets in single_quantile_3 are low
maximum_targets %>% filter(problem=="single_quantile_3")
## Draw and histogram to show PU classification and demand points
problems[[1]]@data@attribute.spaces[[10]]@spaces[[1]]@planning.unit.points@coords %>% as.vector -> pu_coords 
problems[[1]]@data@attribute.spaces[[10]]@spaces[[1]]@demand.points@coords %>% as.vector -> dp_coords
### Re-run Classification 
library(rgeoda)
library(DescTools)
source(paste0(getwd(),"/functions/split.taxon.R"))
g_rast <- rast(paste0(getwd(),"/Results May_2023/Diplodus_allAxes_8068.grd"))
g_rast_values <- terra::extract(g_rast,pus_centroid)
load("Planning_units.RData")
res_split.taxon <- split.taxon(x = g_rast_values[,10],
                               num_classes = 3,
                               class_method = "quantile",
                               rij_taxon = pull(pus,"Diplodus_sargus"),
                               name_feat = paste0("Diplodus_a",sprintf("%02d",10)),
                               return_class_midpoint = T)
### Histogram bars will be coloured according to the class of a PU
res_split.taxon[[2]] %>% apply(1,which.max) -> color_pus
### Create a new data,frame with colors and only the occupied PUs
pus %>% st_drop_geometry %>% transmute(Diplodus_sargus = Diplodus_sargus, color = color_pus) %>% filter(Diplodus_sargus == 1) %>%
    mutate(pu_coords = pu_coords) -> pu_coords_color
### Plot the histogram
png(paste0("Example_dp_raptr_",species,".png"),width=15,height=10,units="cm",res=300)
pu_coords_color %>% filter(color==1) %>% pull(pu_coords) %>%
    hist(breaks=seq(-11.75,15.5,0.25), border=NA,
         xlab="PU coordinates (PCA axis scores) in genetic space #10", main="", ylab="Number of PUs", col="navy")
pu_coords_color %>% filter(color==2) %>% pull(pu_coords) %>% hist(breaks=seq(-11.75,15.5,0.25), border=NA, add=T, col="orange")
pu_coords_color %>% filter(color==3) %>% pull(pu_coords) %>% hist(breaks=seq(-11.75,15.5,0.25), border=NA, add=T, col="darkgreen")
matrix(c(dp_coords,rep(-2.5,3)),nrow=3) %>% points(pch=24,col="black", bg=c("navy","orange","darkgreen"), cex=0.80)
dev.off()
