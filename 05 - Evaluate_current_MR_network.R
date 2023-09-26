# Evaluate how the current MR network satisfy the targets of the problems defined by prioritizr and raptr

rm(list=ls())

species <- "Mullus"

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

