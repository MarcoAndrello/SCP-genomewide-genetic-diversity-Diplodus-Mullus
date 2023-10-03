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

load(paste0("Results_prioritizr_",species,".RData"))

#####################################################################################
# (1) pairwise distance between PU selection frequencies
#####################################################################################
# i) I calculate a selection frequency for each of the two scenarios
# ii) I use the 2 selection frequencies to calculate a Jaccard distance
# iii) I assess significance by permuting solution vectors between the 2 scenarios 
jaccard_distance <- pvalue <- array(NA,dim=c(length(results),length(results)))
dimnames(jaccard_distance) <- dimnames(pvalue) <- list("1" = names(results), "2" = names(results))

par(mar=c(1,1,4,1),mfrow=c(4,4))
for (i.res in 1 : (length(results)-1)) {
    for (j.res in (i.res+1) : length(results)) {
        cat(i.res,j.res,"\n"); flush.console()
        twogroups <- cbind(
            results[[i.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100))),
            results[[j.res]] %>% st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
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
# (started at 9:13)
save(jaccard_distance, pvalue, file=paste0("Jaccard_distance_prioritizr_",species,".RData"))
# Plot a tree showing distances between solutions
load(paste0("Jaccard_distance_prioritizr_",species,".RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png(paste0("hclust_prioritizr_",species,".png"),width=9,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main=paste("prioritizr",species),xlab="",sub="")
dev.off()



#####################################################################################
# (2) Targets met by solutions found with different scenarios
#####################################################################################
# For each scenario, I analyze whether the other scenarios meet its target,by calculating relative targets over conservation features
amount_held <- array(NA,c(length(problems),length(results),100))
i.prob <- j.prob <- 1
for (i.prob in 1 : length(problems)) {
    for (j.prob in 1 : length(results)) {
        for (i.sol in 1 : 100) {
            if (i.sol %% 25 == 0) {
                cat(i.prob,j.prob,i.sol,"\n");
                flush.console()
            }
            problems[[i.prob]] %>%
                eval_target_coverage_summary(select(results[[j.prob]],paste0("solution_",i.sol))) %>%
                filter(relative_target > 0) %>%
                pull(relative_held) %>%
                mean() ->
                amount_held[i.prob,j.prob,i.sol]
        }
    }
}
save(amount_held,file=paste0("Amount_held_prioritizr_",species,".RData"))
load(paste0("Amount_held_prioritizr_",species,".RData"))

## Reshaping amount_held into a dataframe
amount_held_list <- list()
for(i.prob in 1 : 12) {
    amount_held_list[[i.prob]] <-
        amount_held[i.prob,,] %>%
        as.data.frame() %>%
        mutate(solution = factor(names(results),levels=names(results))) %>%
        pivot_longer(cols=1:100) %>%
        mutate(amount_held = value, problem = names(problems)[i.prob]) %>%
        select(-value)
}
amount_held <- bind_rows(amount_held_list)
amount_held$problem <- factor(amount_held$problem,
                            levels=names(problems))
rm(amount_held_list)
# Summarise by problem and scenario
amount_held %>% group_by(solution, problem) %>% summarise(mean_amount_held=mean(amount_held)) %>% 
    mutate(amount_held=mean_amount_held) %>% select(-mean_amount_held) -> mean_amount_held
## Define method and num_clasees column
mean_amount_held$method = c(rep("quantile",12*3),
rep("equal",12*3),
rep("natural",12*3),
rep("sd",12*1),
rep("kmeans",12*2))
mean_amount_held$method = factor(mean_amount_held$method,
                levels=c("quantile","equal","natural","sd","kmeans"))
mean_amount_held$num_classes=c(rep(
    rep(c("low","medium","high"),each=12),
    3),
    rep("medium",12),
    rep("low",12),
    rep("high",12)
)
mean_amount_held$num_classes = factor(mean_amount_held$num_classes,levels=c("low","medium","high"))

# Plot
png(paste0("amount_held_prioritizr_",species,".png"),width=20,height=10,units="cm",res=300)
theme_set(theme_classic())
ggplot(mean_amount_held,aes(x=problem,y=amount_held)) +
    geom_point(aes(colour=method, shape=num_classes),position=position_dodge(width=0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    geom_hline(yintercept=0.15,linetype="dotted") +
    scale_fill_brewer(type="qual") +
    ylim(0.1,0.2005)
dev.off()


# Raster
# png(paste0("amount_held_prioritizr_",species,".png"),width=11,height=7.5,units="cm",res=300)
# ggplot(mean_amount_held,aes(x=problem,y=scenario)) +
#     geom_raster(aes(fill=amount_held)) +
#     scale_fill_gradient(low="red", high="blue", limits=c(0,0.20)) +
#     theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
# dev.off()

# Old: boxplots in panels
# # Reshaping shortfall into a dataframe
# shortfall_list <- list()
# for(i.prob in 1 : 12) {
#     shortfall_list[[i.prob]] <-
#     shortfall[i.prob,,] %>%
#     as.data.frame() %>%
#     mutate(scenario=factor(names(results),levels=names(results))) %>%
#     pivot_longer(cols=1:100) %>%
#     mutate(shortfall = value, problem = names(problems)[i.prob]) %>%
#     select(-value)
# }
# shortfall <- bind_rows(shortfall_list)
# shortfall$problem <- factor(shortfall$problem,
#                             levels=names(problems))
# rm(shortfall_list)
# # Summarise by problem and scenario
# shortfall %>% group_by(scenario, problem) %>% summarise(mean_shortfall=mean(shortfall)) -> mean_shortfall
# # Plot
# png(paste0("Shortfall_prioritizr_",species,".png"),width=11,height=7.5,units="cm",res=300)
# ggplot(mean_shortfall,aes(x=problem,y=scenario)) +
#     geom_raster(aes(fill=mean_shortfall)) +
#     theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
# dev.off()
# # Plot
# png(paste0("Shortfall_prioritizr_",species,"_distribution.png"),width=15,height=20,units="cm",res=300)    
# theme_set(theme_classic())
# ggplot(shortfall,aes(x=scenario,y=shortfall)) +
#     geom_boxplot() +
#     geom_hline(yintercept=2,linetype="dotted") +
#     facet_wrap(vars(problem),ncol=3) +
#     theme(axis.text.x = element_text(angle = 90))
# dev.off()


#####################################################################################
# (3) Conservation cost (number of selected PUs)
#####################################################################################
cost <- array(NA,c(length(problems),100))
i.prob <- 1
for (i.prob in 1 : length(problems)) {
    for (i.sol in 1 : 100) {
        if (i.sol %% 25 == 0) {
            cat(i.prob,i.sol,"\n");
            flush.console()
        }
        problems[[i.prob]] %>%
            # eval_cost_summary(select(results[[i.prob]],paste0("solution_",i.sol))) %>%
            # pull(cost) ->
            eval_n_summary(select(results[[i.prob]],paste0("solution_",i.sol))) %>% 
            pull(n) ->
            cost[i.prob,i.sol]
    }
}
save(cost,file=paste0("Cost_prioritizr_",species,".RData"))
load(paste0("Cost_prioritizr_",species,".RData"))
# Reshaping cost into a dataframe
cost %>% t() -> cost
colnames(cost) <- names(problems)
cost <- as.data.frame(cost)
cost %>%
    pivot_longer(col=names(problems),
                 names_to="problem",
                 values_to="cost",
                 cols_vary="slowest") ->
    cost
cost$problem <- factor(cost$problem,
                            levels=names(problems))


## Define method and num_clasees column
cost$method = c(rep("quantile",100*3),
                            rep("equal",100*3),
                            rep("natural",100*3),
                            rep("sd",100*1),
                            rep("kmeans",100*2))
cost$method = factor(cost$method,
                                 levels=c("quantile","equal","natural","sd","kmeans"))
cost$num_classes=c(rep(
    rep(c("low","medium","high"),each=100),
    3),
    rep("medium",100),
    rep("low",100),
    rep("high",100)
)
cost$num_classes = factor(cost$num_classes,levels=c("low","medium","high"))
# Cost required to protect the species distribution only
cost_species_distribution <- ifelse(species == "Diplodus", 340, 544)
# Plot
png(paste0("Cost_prioritizr_",species,".png"),width=7.5,height=7.5,units="cm",res=300)    
theme_set(theme_classic())
ggplot(cost,aes(x=problem,y=cost)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
    geom_hline(yintercept=cost_species_distribution,linetype="dotted") +
    ggtitle(paste0(paste("prioritizr",species)))
dev.off()

mod <- lm(cost~method*num_classes,data=cost)


#####################################################################################
# Other code
#####################################################################################

######
# Compare distances
par(mfrow=c(2,2))
vegdist_method <- c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower",
                    "horn", "mountford",  "binomial", "chord", "hellinger",  "robust.aitchison") # "aitchison","chisq", "chao", "cao", "mahalanobis", "morisita", "raup",
for (i.method in 1 : length(vegdist_method)){
    selection_frequency_1 <- rowSums(twogroups[,1:100]) / 100
    selection_frequency_2 <- rowSums(twogroups[,101:200]) / 100
    vegdist(rbind(selection_frequency_1,selection_frequency_2),method=vegdist_method[i.method]) %>% as.numeric() -> distance_observed
    cat(vegdist_method[i.method],distance_observed,"\n"); flush.console()
    # Permute columns
    num_perm <- 100
    distance_permuted <- rep(NA,num_perm)
    for (i.perm in 1 : num_perm) {
        if (i.perm %% 100 == 0) {
            cat(i.perm,"of",num_perm,"\n")
            flush.console()
        }
        twogroups_permuted <- twogroups[,sample(200)]
        selection_frequency_1 <- rowSums(twogroups_permuted[,1:100]) / 100
        selection_frequency_2 <- rowSums(twogroups_permuted[,101:200]) / 100
        vegdist(rbind(selection_frequency_1,selection_frequency_2),method=vegdist_method[i.method]) %>% as.numeric() -> distance_permuted[i.perm]
    }
    hist(distance_permuted,main=vegdist_method[i.method],
         xlim=c(min(distance_permuted),
                max((distance_observed+1), max(distance_permuted))))
    abline(v=distance_observed)
}

