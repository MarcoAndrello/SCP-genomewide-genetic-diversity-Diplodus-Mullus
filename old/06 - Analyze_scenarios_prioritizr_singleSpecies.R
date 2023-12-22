# Assess scenarios prioritizr
# (1) pairwise distance between PU selection frequencies
# (2) Targets met by solutions found with different scenarios
# (3) Conservation cost

rm(list=ls())

species <- "Mullus"

library(tidyverse)
library(sf)
library(ade4)
library(vegan)
library(prioritizr)
library(raptr)
# library(gridExtra)

load(paste0(getwd(),"/Results/Results_prioritizr_",species,".RData"))
load(paste0(getwd(),"/Results/Results_raptr_",species,".RData")) # res_gs


#####################################################################################
# (1) pairwise distance between PU selection frequencies
#####################################################################################
# i) I calculate a selection frequency for each of the two scenarios
# ii) I use the 2 selection frequencies to calculate a Jaccard distance
# iii) I assess significance by permuting solution vectors between the 2 scenarios 

jaccard_distance <- pvalue <- array(NA,dim=c(length(results)+1,length(results)+1))
dimnames(jaccard_distance) <- dimnames(pvalue) <- list("1" = c(names(results),"raptr"), "2" = c(names(results),"raptr"))

par(mar=c(1,1,4,1),mfrow=c(4,4))
# i.res <- 1; j.res <- 13
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
        selection_frequency_2 <- rowSums(twogroups[,101:200]) / 100
        vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_observed
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
        jaccard_distance[i.res,j.res] <- distance_observed
        pvalue[i.res,j.res] <- 1 - (ecdf(distance_permuted)(distance_observed))
    }
}

save(jaccard_distance, pvalue, file=paste0(getwd(),"/Results/Jaccard_distance_",species,".RData"))
# Plot a tree showing distances between solutions
load(paste0(getwd(),"/Results/Jaccard_distance_",species,".RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png(paste0("Figures/hclust_prioritizr_",species,".png"),width=9,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main=paste("prioritizr",species),xlab="",sub="")
dev.off()



#####################################################################################
# (2) Space held by solutions found with different scenarios
#####################################################################################
list_space_held <- list()
# Solutions found with prioritizr
for (i.prob in 1 : length(results)) { 
    space_held_solution <- vector()
    for (i.sol in 1 : 100) {
        if (i.sol %% 25 == 0) {
            cat(i.prob,i.sol,"\n");
            flush.console()
        }
        selections <- which(results[[i.prob]] %>% pull(paste0("solution_",i.sol))==1)
        res_updated <- update(res_gs, b = selections)
        space_held_solution[i.sol] <- space.held(res_updated, y=1, space=NULL) %>% as.vector()
    }
    data.frame(space_held = space_held_solution,
               solution = names(problems)[[i.prob]],
               sol = as.integer(1:100)) ->
        list_space_held[[i.prob]]
}

# Solutions found with raptr
space_held_solution <- space.held(res_gs, y=NULL,space=NULL) %>% as.vector()
data.frame(space_held = space_held_solution,
           solution = "raptr",
           sol = as.integer(1:100)) ->
    list_space_held[[13]]
    
space_held <- bind_rows(list_space_held)
space_held$solution <- factor(space_held$solution, levels=c(names(results),"raptr"))
save(space_held, file="space_held_mullus.RData")

# Plot
theme_set(theme_classic())
png(paste0("Figures/Space_held_",species,".png"),width=15,height=10,units="cm",res=300)
ggplot(space_held) +
    geom_boxplot(aes(x=solution, y=space_held)) +
    # geom_violin(aes(x=solution, y=space_held), fill="black") #+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    ylim(0.7,1) +
    geom_hline(yintercept=0.75, linetype = "dotted") +
    ggtitle(species)
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
cost_solution <- res_gs@results@summary$Cost
data.frame(cost = cost_solution,
           solution = "raptr",
           sol = as.integer(1:100)) ->
    list_cost[[13]]

cost <- bind_rows(list_cost)
cost$solution <- factor(cost$solution, levels=c(names(results),"raptr"))

# Number of cells and cost required to protect the species distribution only
n_species_distribution <- ifelse(species == "Diplodus", 338, 542) # before was 340, 544
#cost_species_distribution ### fare scenario con prioritizr


# Plot
theme_set(theme_classic())
png(paste0("Figures/Cost_",species,".png"),width=15,height=10,units="cm",res=300)
ggplot(cost) +
    geom_boxplot(aes(x=solution, y=cost)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    # geom_hline(yintercept=cost_species_distribution,linetype="dotted") +
    ggtitle(species) + ylab("Cost (k€)")
dev.off()


## Significant differences?
## Define method and num_classes column
cost$method = c(rep("quantile",100*3),
                rep("equal",100*3),
                rep("natural",100*3),
                rep("sd",100*1),
                rep("kmeans",100*2),
                rep("raptr",100*1)
)
cost$method = factor(cost$method,
                                 levels=c("quantile","equal","natural","sd","kmeans","raptr"))
cost$num_classes=c(rep(
    rep(c("low","medium","high"),each=100),
    3),
    rep("medium",100),
    rep("low",100),
    rep("high",100),
    rep("raptr-gs",100)
)
cost$num_classes = factor(cost$num_classes,levels=c("low","medium","high","raptr-gs"))
mod <- lm(cost~method*num_classes,data=cost)
