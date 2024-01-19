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
library(tmap)
library(rnaturalearth)
library(RColorBrewer)

load("Results/Problem_definition_raptr.RData")
load("Results/Results_prioritizr.RData")
load("Results/Results_raptr_50gs.RData")
load("Results/Results_raptr_20gs.RData")


#####################################################################################
# (1) pairwise distance between PU selection frequencies
#####################################################################################
# i) I calculate a selection frequency for each of the two scenarios
# ii) I use the 2 selection frequencies to calculate a Jaccard distance
# iii) I assess significance by permuting solution vectors between the 2 scenarios 

jaccard_distance <- pvalue <- array(NA,dim=c(length(results_prioritizr)+2,length(results_prioritizr)+2))
dimnames(jaccard_distance) <- dimnames(pvalue) <- list("1" = c(names(results_prioritizr),"raptr_50perc","raptr_20perc"),
                                                       "2" = c(names(results_prioritizr),"raptr_50perc","raptr_20perc"))

for (i.res in 1 : (nrow(jaccard_distance)-1)) {
    for (j.res in (i.res+1) : nrow(jaccard_distance) ) {
        cat(i.res,j.res,"\n"); flush.console()
        if (i.res < 14) {
            first_group <- results_prioritizr[[i.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
        }
        if(i.res == 14) {
            first_group <- res_50gs[[1]]@results@selections
            for (i in 2 : 20) first_group <- rbind(first_group,
                                                   res_50gs[[i]]@results@selections)
            first_group <- t(first_group)
        }
        if(i.res == 15) {
            first_group <- res_20gs[[1]]@results@selections
            for (i in 2 : 20) first_group <- rbind(first_group,
                                                   res_20gs[[i]]@results@selections)
            first_group <- t(first_group)
         
        }
        if (j.res < 14) {
            second_group <- results_prioritizr[[j.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
        }
        if(j.res == 14) {
            second_group <- res_50gs[[1]]@results@selections
            for (i in 2 : 20) second_group <- rbind(second_group,
                                                   res_50gs[[i]]@results@selections)
            second_group <- t(second_group)

        }
        if(j.res == 15) {
            second_group <- res_20gs[[1]]@results@selections
            for (i in 2 : 20) second_group <- rbind(second_group,
                                                   res_20gs[[i]]@results@selections)
            second_group <- t(second_group)

        }

        twogroups <- cbind(first_group, second_group)

        selection_frequency_1 <- rowSums(first_group) / ncol(first_group)
        selection_frequency_2 <- rowSums(second_group) / ncol(second_group)
        vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_observed
        jaccard_distance[i.res,j.res] <- distance_observed

        # Permute columns
        num_perm <- 1000
        distance_permuted <- rep(NA,num_perm)
        for (i.perm in 1 : num_perm) {
            if (i.perm %% 100 == 0) {
                cat(i.perm,"of",num_perm,"\n")
                flush.console()
            }
            twogroups_permuted <- twogroups[,sample(ncol(twogroups))]
            selection_frequency_1 <- rowSums(twogroups_permuted[,1:ncol(first_group)]) / ncol(first_group)
            selection_frequency_2 <- rowSums(twogroups_permuted[,(ncol(first_group)+1):(ncol(first_group)+ncol(second_group))]) / ncol(second_group)
            vegdist(rbind(selection_frequency_1,selection_frequency_2),method="jaccard") %>% as.numeric() -> distance_permuted[i.perm]
        }
        pvalue[i.res,j.res] <- 1 - (ecdf(distance_permuted)(distance_observed))
    }
}

jaccard_distance %>% as.vector %>% summary
save(jaccard_distance, pvalue, file="Results/Jaccard_distance.RData")

# Plot a tree showing distances between solutions
load(paste0(getwd(),"/Results/Jaccard_distance.RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png("Figures/hclust.png",width=15,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main="Distance between solutions",xlab="",sub="")
dev.off()


#####################################################################################
# (2) Space held by solutions found with different scenarios
#####################################################################################
list_space_held <- list()
# Solutions found with prioritizr
# Beware: takes about 2 hrs per i.prob
for (i.prob in 1 : length(results_prioritizr)) { 
    space_held_Diplodus <- space_held_Mullus <- vector()
    for (i.sol in 1 : 100) {
        cat(i.prob,i.sol,"\n")
        flush.console()
        selections <- which(results_prioritizr[[i.prob]] %>% pull(paste0("solution_",i.sol))==1)
        res_updated <- update(prob_gs, b = selections)
        space_held_Diplodus[i.sol] <- space.held(res_updated, y=1, species=1) %>% as.vector()
        space_held_Mullus[i.sol] <- space.held(res_updated, y=1, species=2) %>% as.vector()
    }
    data.frame(space_held_Diplodus = space_held_Diplodus,
               space_held_Mullus = space_held_Mullus,
               solution = names(results_prioritizr)[[i.prob]],
               sol = as.integer(1:100)) ->
        list_space_held[[i.prob]]
}
save(list_space_held, file=paste0("list_space_held_12e13.RData"))

# # Solutions found with raptr
# res_50gs
space_held_Diplodus <- space_held_Mullus <- vector()
for (i.sol in 1 : 20) {
    cat(i.sol,"\n")
    flush.console()
    res_50gs[[i.sol]]@results@selections %>%
        t %>%
        as.vector %>%
        `==`(1) %>%
        which ->
        selections
    update(prob_gs, b = selections) ->
        res_updated
    space.held(res_updated, y=1, species=1) %>% as.vector() ->
        space_held_Diplodus[i.sol]
    space.held(res_updated, y=1, species=2) %>% as.vector() ->
        space_held_Mullus[i.sol]
}
data.frame(space_held_Diplodus = space_held_Diplodus,
           space_held_Mullus = space_held_Mullus,
           solution = "raptr_50perc",
           sol = as.integer(1:20)) ->
list_space_held[[14]]


# res_20gs
space_held_Diplodus <- space_held_Mullus <- vector()
for (i.sol in 1 : 20) {
    cat(i.sol,"\n")
    flush.console()
    res_20gs[[i.sol]]@results@selections %>%
        t %>%
        as.vector %>%
        `==`(1) %>%
        which ->
        selections
    update(prob_gs, b = selections) ->
        res_updated
    space.held(res_updated, y=1, species=1) %>% as.vector() ->
        space_held_Diplodus[i.sol]
    space.held(res_updated, y=1, species=2) %>% as.vector() ->
        space_held_Mullus[i.sol]
}
data.frame(space_held_Diplodus = space_held_Diplodus,
           space_held_Mullus = space_held_Mullus,
           solution = "raptr_20perc",
           sol = as.integer(1:20)) ->
    list_space_held[[15]]

save(list_space_held,file="Results/List_space_held.RData")
load("Results/List_space_held.RData")
space_held <- bind_rows(list_space_held)
space_held$solution <- factor(space_held$solution, levels=unique(space_held$solution))

# Plot
space_held %>% rename(Diplodus = space_held_Diplodus, Mullus = space_held_Mullus) %>% 
    pivot_longer(cols=c(Diplodus, Mullus)) %>%
    rename(species=name, space_held=value) -> space_held_longer
theme_set(theme_classic())
png(paste0("Figures/Space_held.png"),width=15,height=10,units="cm",res=300)
ggplot(space_held_longer) +
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
for (i.prob in 1 : length(problems_prioritizr) ) {
    cost_solution <- vector()
    for (i.sol in 1 : 100) {
        if (i.sol %% 25 == 0) {
            cat(i.prob,i.sol,"\n");
            flush.console()
        }
        problems_prioritizr[[i.prob]] %>%
            eval_cost_summary(select(results_prioritizr[[i.prob]],paste0("solution_",i.sol))) %>% pull(cost) ->
            cost_solution[i.sol]
    }
    data.frame(cost = cost_solution,
               solution = names(problems_prioritizr)[[i.prob]],
               sol = as.integer(1:100)) ->
        list_cost[[i.prob]]
}

# Cost of raptr solutions
# res_50gs
cost_solution <- vector()
for (i.sol in 1 : 20) {
    res_50gs[[i.sol]]@results@summary$Cost ->
        cost_solution[i.sol]
}
data.frame(cost = cost_solution,
           solution = "raptr_50perc",
           sol = as.integer(1:20)) ->
    list_cost[[14]]
# res_20gs
cost_solution <- vector()
for (i.sol in 1 : 20) {
    res_20gs[[i.sol]]@results@summary$Cost ->
        cost_solution[i.sol]
}
data.frame(cost = cost_solution,
           solution = "raptr_20perc",
           sol = as.integer(1:20)) ->
    list_cost[[15]]

cost <- bind_rows(list_cost)
cost$solution <- factor(cost$solution, levels=unique(cost$solution))
cost$cost <- cost$cost / 1000 # transforms kilo-euro into Million-euros
save(cost,file="Results/Cost.RData")

cost %>% group_by(solution) %>% summarise(mean_cost=mean(cost)) %>% arrange(mean_cost)

# Plot
theme_set(theme_classic())
png(paste0("Figures/Cost.png"),width=15,height=10,units="cm",res=300)
ggplot(cost) +
    geom_boxplot(aes(x=solution, y=cost)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    geom_hline(yintercept=1047,linetype="dotted") +
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

# Relationships between space_held and conservation cost
space_held %>% mutate(cost= cost$cost) %>% group_by(solution) %>%
    summarise(mean_cost=mean(cost), mean_space_held=mean(space_held_Diplodus)) -> data
ggplot(data,aes(x=mean_space_held, y=mean_cost, col=solution)) +
    geom_point()


#####################################################################################
# (4) Map of selection frequencies
#####################################################################################
# Load planning units and countries shapefiles
load("Planning_units.RData")
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(pus)) -> countries

# Loop on planning problems to plot one map per problem
for (i.res in 14 : 15) {
    cat(i.res,"\n"); flush.console()
    if (i.res < 14) {
        selections <- results_prioritizr[[i.res]] %>%
            st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
        names(problems_prioritizr)[i.res] -> main.title
    }
    if(i.res == 14) {
        selections <- res_50gs[[1]]@results@selections
        for (i in 2 : 20) selections <- rbind(selections,
                                               res_50gs[[i]]@results@selections)
        selections <- t(selections)
        main.title <- "raptr_50perc"
        selections -> best_selections
    }
    if(i.res == 15) {
        selections <- res_20gs[[1]]@results@selections
        for (i in 2 : 20) selections <- rbind(selections,
                                              res_20gs[[i]]@results@selections)
        selections <- t(selections)
        # selections <- rbind(results_raptr[[6]]@results@selections,
        #                     results_raptr[[7]]@results@selections,
        #                     results_raptr[[8]]@results@selections,
        #                     results_raptr[[9]]@results@selections,
        #                     results_raptr[[10]]@results@selections) %>% t()
        main.title <- "raptr_20perc"
    }
    selection_frequency <- rowSums(selections) / ncol(selections)
    
# Add selection frequency 
    pus %>%
        mutate(sel_frequency =
                   # Divide selection frequency into breaks
                   selection_frequency %>%
                   cut(breaks=seq(0,1,0.2),include.lowest=T,right=F) %>%
                   # Make factor with levels
                   factor(levels = c("Existing","[0,0.2)","[0.2,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1]"))
        ) %>%
        # Replace selection_frequency of existing reserves with level "Existing"
        mutate(selection_frequency = replace(sel_frequency,
                                             which(pus$status_0 == 1), # These are the existing reserves
                                             "Existing")
        ) -> pus_map
    
    png(paste0("Figures/Maps_selections/Map_selection_frequency_prob",i.res,".png"),
        width=20,height=8,res=600,units="cm")
    print(
        tm_shape(pus_map) +
        tm_fill("selection_frequency", palette=c("gold",brewer.pal(5,"Greens")),
                legend.is.portrait = F, legend.show = F) +
        tm_shape(countries, bbox = res) +
        tm_polygons(col="lightgray", lwd=0.2) +
        tm_layout(main.title=main.title, main.title.position = "center")
    )
    dev.off()
    
    # Plot map of the best solutions (that of raptr_50)
    if (i.res == 14) {
        png("Figures/Map_selection_frequency.png",
            width=20,height=15,res=600,units="cm")
        print(
            tm_shape(pus_map) +
                tm_fill("selection_frequency", palette=c("gold",brewer.pal(5,"Greens")),
                        legend.is.portrait = F, legend.show = T) +
                tm_shape(countries, bbox = res) +
                tm_polygons(col="lightgray", lwd=0.2) +
                tm_legend(legend.outside = T, legend.outside.position = "bottom") +
                tm_layout(main.title=main.title, main.title.position = "center")
        )
        dev.off()
        }
    
}

# Plot the legend only
png("Figures/Maps_selections/Map_selection_legend.png",
    width=20,height=8,res=600,units="cm")
print(
    tm_shape(pus_map) +
        tm_fill("selection_frequency", palette=c("gold",brewer.pal(5,"Greens")),
                legend.is.portrait = F, legend.show = T) +
        tm_legend(legend.outside = T, legend.outside.position = "bottom") +
        tm_layout(legend.only = T)
)
dev.off()







