# Evaluate how the current MR network satisfy the targets of the problems defined by prioritizr and raptr

rm(list=ls())

species <- "Diplodus"

library(tidyverse)
library(sf)
library(prioritizr)
library(raptr)

# Load planning units
load("Planning_units.RData")


#####################################################################################
# Prioritizr, targets met 
#####################################################################################
load(paste0("Results_prioritizr_",species,".RData"))
rm(results)
shortfall <- vector()
i.prob <- 1
for (i.prob in 1 : length(problems)) {
    eval_target_coverage_summary(problems[[i.prob]],
                                 select(pus, status_0)) %>%
        pull(relative_shortfall) %>%
        mean ->
        shortfall[[i.prob]]
}    
shortfall <- tibble(problem=names(problems),shortfall=shortfall)
shortfall$problem <- factor(shortfall$problem,
                            levels=names(problems))
png(paste0("Evaluate_MR_prioritizr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)
theme_set(theme_classic())
ggplot(shortfall,aes(x=problem,y=shortfall)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    ylim(0,0.15)
dev.off()
rm(problems, shortfall, i.prob)
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
    distance_to_maximum[i.prob] <- mean(maximum_targets_problem - space_held_feature)
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

distance_to_maximum <- tibble(problem=names(problems),distance_to_maximum=distance_to_maximum)
distance_to_maximum$problem <- factor(distance_to_maximum$problem,
                              levels=names(problems))
png(paste0("Evaluate_MR_raptr_distance_",species,".png"),width=7.5,height=7.5,units="cm",res=300)
theme_set(theme_classic())
ggplot(distance_to_maximum ,aes(x=problem,y=distance_to_maximum )) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) 
dev.off()
rm(problems, results, space_held, distance_to_maximum, i.prob)






# ## OLD CODE PRIORITIZR
# # For each scenario, I analyze whether the current MR network meets its target, by calculating total shortfall over conservation features
# shortfall_list <- list()
# i.prob <- 1
# for (i.prob in 1 : length(problems)) {
#     eval_target_coverage_summary(problems[[i.prob]],
#                                  select(pus, status_0)) -> 
#         shortfall_list[[i.prob]]
# }
# names(shortfall_list) <- names(problems)
# bind_rows(shortfall_list, .id="problem") %>%
#     select(problem, feature, relative_shortfall) -> shortfall
# shortfall$problem <- factor(shortfall$problem,
#                             levels=names(problems))
# rm(shortfall_list)
# # Plot
# png(paste0("Evaluate_MR_prioritizr_",species,".png"),width=20,height=10,units="cm",res=300)    
# theme_set(theme_classic())
# ggplot(shortfall,aes(x=problem, y=relative_shortfall)) +
#     geom_dotplot(binaxis="y",stackdir="center", binwidth=0.001) +
#     theme(axis.text.x = element_text(angle = 90))
# dev.off()
# rm(problems, shortfall, i.prob)