# Assess scenarios prioritizr
# (1) pairwise distance between PU selection frequencies
# (2) Targets met by solutions found with different scenarios
# (3) Conservation cost

rm(list=ls())

library(tidyverse)
library(sf)
library(ade4)
library(vegan)
library(prioritizr)
library(raptr)

load(paste0(getwd(),"/Results/Results_prioritizr.RData"))
load(paste0(getwd(),"/Results/Results_raptr.RData"))


#####################################################################################
# (1) pairwise distance between PU selection frequencies
#####################################################################################
# i) I calculate a selection frequency for each of the two scenarios
# ii) I use the 2 selection frequencies to calculate a Jaccard distance
# iii) I assess significance by permuting solution vectors between the 2 scenarios 

jaccard_distance <- pvalue <- array(NA,dim=c(length(results)+1,length(results)+1))
dimnames(jaccard_distance) <- dimnames(pvalue) <- list("1" = c(names(results),"raptr"), "2" = c(names(results),"raptr"))

par(mar=c(1,1,4,1),mfrow=c(4,4))
i.res <- 1; j.res <- 13
for (i.res in 1 : length(results)) {
    for (j.res in (i.res+1) : (length(results)+1) ) {
        cat(i.res,j.res,"\n"); flush.console()
        twogroups <- cbind(
            results[[i.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100))),
            if (j.res < 13) {
                results[[j.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
            } else {
                res_gs@results@selections %>% t()
            }
        )
        selection_frequency_1 <- rowSums(twogroups[,1:100]) / 100
        selection_frequency_2 <- if(j.res < 13) rowSums(twogroups[,101:200])/100 else twogroups[,101]
        vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_observed
        jaccard_distance[i.res,j.res] <- distance_observed
        
        if(j.res == 13) next
        # Permute columns
        num_perm <- 1000
        distance_permuted <- rep(NA,num_perm)
        for (i.perm in 1 : num_perm) {
            if (i.perm %% 100 == 0) {
                cat(i.perm,"of",num_perm,"\n")
                flush.console()
            }
            twogroups_permuted <- twogroups[,sample(200)]
            selection_frequency_1 <- rowSums(twogroups_permuted[,1:100]) / 100
            selection_frequency_2 <- rowSums(twogroups_permuted[,101:200]) / 100
            vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_permuted[i.perm]
        }
        hist(distance_permuted,
             main = names(results)[c(i.res,j.res)],
             xlim=c(min(distance_permuted),
                    max((distance_observed+1), max(distance_permuted))))
        abline(v=distance_observed,col="red")
        
        pvalue[i.res,j.res] <- 1 - (ecdf(distance_permuted)(distance_observed))
    }
}
save(jaccard_distance, pvalue, file=paste0(getwd(),"/Results/Jaccard_distance.RData"))

# Plot a tree showing distances between solutions
load(paste0(getwd(),"/Results/Jaccard_distance.RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png(paste0("Figures/hclust.png"),width=11,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main="Distance between solutions",xlab="",sub="")
dev.off()



#####################################################################################
# (2) Space held by solutions found with different scenarios
#####################################################################################
list_space_held <- list()
# Solutions found with prioritizr
for (i.prob in 1 : length(results)) { 
    space_held_Diplodus <- space_held_Mullus <- vector()
    for (i.sol in 1 : 100) {
        if (i.sol %% 25 == 0) {
            cat(i.prob,i.sol,"\n");
            flush.console()
        }
        selections <- which(results[[i.prob]] %>% pull(paste0("solution_",i.sol))==1)
        res_updated <- update(res_gs, b = selections)
        space_held_Diplodus[i.sol] <- space.held(res_updated, y=1, species=1) %>% as.vector()
        space_held_Mullus[i.sol] <- space.held(res_updated, y=1, species=2) %>% as.vector()
    }
    data.frame(space_held_Diplodus = space_held_Diplodus,
               space_held_Mullus = space_held_Mullus,
               solution = names(problems)[[i.prob]],
               sol = as.integer(1:100)) ->
        list_space_held[[i.prob]]
}


# # Solutions found with raptr
# space_held_Diplodus <- space.held(res_gs, y=1, species=1) %>% as.vector()
# space_held_Mullus <- space.held(res_gs, y=1, species=1) %>% as.vector()
# data.frame(space_held = space_held_solution,
#            solution = "raptr",
#            sol = as.integer(1:100)) ->
#     list_space_held[[13]]
    
space_held <- bind_rows(list_space_held)
space_held$solution <- factor(space_held$solution, levels=names(results))
save(space_held, file="Results/Space_held.RData")

# Plot
load("Results/Space_held.RData")
space_held %>% rename(Diplodus = space_held_Diplodus, Mullus = space_held_Mullus) %>% 
    pivot_longer(cols=c(Diplodus, Mullus)) %>%
    rename(species=name, space_held=value) -> space_held
theme_set(theme_classic())
png(paste0("Figures/Space_held.png"),width=15,height=10,units="cm",res=300)
ggplot(space_held) +
    geom_boxplot(aes(x=solution, y=space_held, fill=species, colour=species)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    ylim(0.85,1) +
    geom_hline(yintercept=0.95, linetype = "dotted")
dev.off()


#####################################################################################
# (3) Conservation cost
#####################################################################################
list_cost <- list()
# Cost of prioritizr solutions
i.prob <- 1
for (i.prob in 1 : length(problems) ) {
    cost_solution <- vector()
    for (i.sol in 1 : 100) {
        if (i.sol %% 25 == 0) {
            cat(i.prob,i.sol,"\n");
            flush.console()
        }
        problems[[i.prob]] %>%
            eval_cost_summary(select(results[[i.prob]],paste0("solution_",i.sol))) %>% pull(cost) ->
            cost_solution[i.sol]
    }
    data.frame(cost = cost_solution,
               solution = names(problems)[[i.prob]],
               sol = as.integer(1:100)) ->
        list_cost[[i.prob]]
}
# Cost of raptr solutions
# cost_solution <- res_gs@results@summary$Cost
# data.frame(cost = cost_solution,
#            solution = "raptr",
#            sol = as.integer(1:100)) ->
#     list_cost[[13]]

cost <- bind_rows(list_cost)
cost$solution <- factor(cost$solution, levels=names(results))
cost$cost <- cost$cost / 1000 # transforms kilo-euro into Million-euros

# Number of cells and cost required to protect the species distribution only
# n_species_distribution <- ifelse(species == "Diplodus", 338, 542) # before was 340, 544


# Plot
theme_set(theme_classic())
png(paste0("Figures/Cost.png"),width=15,height=10,units="cm",res=300)
ggplot(cost) +
    geom_boxplot(aes(x=solution, y=cost)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    geom_hline(yintercept=1046,linetype="dotted") +
    geom_hline(yintercept=982,linetype="dashed") +
    ylab("Cost (M€)")
dev.off()


# ## Significant differences?
# ## Define method and num_classes column
# cost$method = c(rep("quantile",100*3),
#                 rep("equal",100*3),
#                 rep("natural",100*3),
#                 rep("sd",100*1),
#                 rep("kmeans",100*2),
#                 rep("raptr",100*1)
# )
# cost$method = factor(cost$method,
#                                  levels=c("quantile","equal","natural","sd","kmeans","raptr"))
# cost$num_classes=c(rep(
#     rep(c("low","medium","high"),each=100),
#     3),
#     rep("medium",100),
#     rep("low",100),
#     rep("high",100),
#     rep("raptr-gs",100)
# )
# cost$num_classes = factor(cost$num_classes,levels=c("low","medium","high","raptr-gs"))
# mod <- lm(cost~method*num_classes,data=cost)
