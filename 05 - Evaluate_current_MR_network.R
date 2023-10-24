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
load(paste0(getwd(),"/Results/Results_prioritizr_",species,".RData"))
rm(results)

# For each scenario, I analyze whether the other scenarios meet its target,by calculating relative targets over conservation features
list_amount_held <- list()
i.prob <- 1
for (i.prob in 1 : length(problems)) {
    eval_target_coverage_summary(problems[[i.prob]],
                                 select(pus, status_0)) %>%
        filter(relative_target > 0) %>%
        select(feature, relative_held) %>%
        mutate(problem = names(problems)[i.prob]) ->
    list_amount_held[[i.prob]]
}
amount_held <- bind_rows(list_amount_held)
amount_held$problem <- factor(amount_held$problem,
                              levels=names(problems))

png(paste0("Evaluate_MR_prioritizr_",species,".png"),width=15,height=10,units="cm",res=300)
theme_set(theme_classic())
ggplot(amount_held,aes(x=problem,y=amount_held)) +
    geom_violin(aes(x=problem, y=relative_held), fill="black") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    geom_hline(yintercept=0.15, linetype="dotted") +
    ggtitle("Existing marine reserves") +
    xlab(paste(species,"planning problems (prioritizr)")) +
    ylim(0,1)
dev.off()
rm(problems, amount_held, i.prob)
#    
    
# Species occurrence held
prob_species_occurrence <- problem(pus,
                                   features = c("Diplodus_sargus","Mullus_surmuletus"),
                                   cost_column = "cost") %>%
    add_relative_targets(0.15) %>%
    add_binary_decisions()
eval_target_coverage_summary(prob_species_occurrence,
                             select(pus, status_0))


#####################################################################################
# Raptr
#####################################################################################
load(paste0(getwd(),"/Results/Results_raptr_",species,".RData"))

selections <- which(pus$status_0 == 1)
list_space_held <- list()
i.prob <- 1
for (i.prob in 1 : length(results)){
    cat(i.prob,"\n"); flush.console()
    # maximum_targets_problem <- maximum.targets(results[[i.prob]])$proportion
    res_updated <- update(results[[i.prob]], b = selections)
    space_held_feature <- space.held(res_updated, y=1, space=NULL) %>%
        t %>%
        as_tibble %>%
        mutate(space_name=row.names(.), problem = names(problems)[i.prob])
    names(space_held_feature)[1] <- "space_held"
    list_space_held[[i.prob]] <- space_held_feature
    if(i.prob > 10) {
        list_space_held[[i.prob]] <- rbind(space_held_feature,space_held_feature,space_held_feature) ### METODO RAPIDO E SPORCO PER CREARE UN VIOLIN PLOT ANCHE QUANDO HO UNA SOLA OSSERVAZIONE
        list_space_held[[i.prob]]$space_held <- list_space_held[[i.prob]]$space_held + rnorm(3,0,0.01)
    }
}
rm(space_held_feature, selections, res_updated)
space_held <- bind_rows(list_space_held)
space_held$problem <- factor(space_held$problem,
                            levels=names(problems))

# space_held %>% group_by(problem) %>% summarise(mean_space_held=mean(space_held)) -> mean_space_held
png(paste0("Evaluate_MR_raptr_",species,".png"),width=15,height=10,units="cm",res=300)
theme_set(theme_classic())
ggplot(space_held,aes(x=problem,y=space_held)) +
    geom_violin(aes(x=problem, y=space_held), fill="black") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    geom_hline(yintercept=0.75, linetype="dotted") +
    ggtitle("Existing marine reserves") +
    xlab(paste(species,"planning problems (raptr)")) +
    ylim(0,1) 
dev.off()


rm(problems, results, space_held, distance_to_maximum, i.prob)

