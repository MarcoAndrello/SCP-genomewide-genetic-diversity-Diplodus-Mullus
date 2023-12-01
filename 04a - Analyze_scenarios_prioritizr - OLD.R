# Assess scenarios prioritizr
# (1) pairwise distance between PU selection frequencies
# (2) Targets met by solutions found with different scenarios
# (3) Conservation cost

rm(list=ls())

species <- "Diplodus"

library(tidyverse)
library(sf)
library(ade4)
library(vegan)
library(prioritizr)
library(gridExtra)

load(paste0(getwd(),"/Results/Results_prioritizr_",species,".RData"))

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

save(jaccard_distance, pvalue, file=paste0(getwd(),"/Results/Jaccard_distance_prioritizr_",species,".RData"))
# Plot a tree showing distances between solutions
load(paste0(getwd(),"/Results/Jaccard_distance_prioritizr_",species,".RData"))
hc <- hclust(as.dist(t(jaccard_distance)))
png(paste0("hclust_prioritizr_",species,".png"),width=9,height=15,units="cm",res=300)
par(mar=c(1,4,2,1))
plot(hc,main=paste("prioritizr",species),xlab="",sub="")
dev.off()



#####################################################################################
# (2) Amount held by solutions found with different scenarios
#####################################################################################

# For each solution, I calculate the amount held in other planning problems
list_amount_held <- list()
plots <- list()
theme_set(theme_classic())
i.prob <- j.prob <- 1
i.prob <- 2
j.prob <- 6
for (i.prob in 1 : length(results)) {
    k <- 1
    for (j.prob in 1 : length(problems)) {
        for (i.sol in 1 : 100) {
            if (i.sol %% 25 == 0) {
                cat("solutions",i.prob,"for problem",j.prob,i.sol,"\n");
                flush.console()
            }
            problems[[j.prob]] %>%
                eval_target_coverage_summary(select(results[[i.prob]],paste0("solution_",i.sol))) %>%
                filter(relative_target > 0) %>%
                select(feature, relative_held) %>%
                mutate(problem = names(problems)[j.prob], solution = names(results)[i.prob], sol = as.integer(i.sol)) ->
                list_amount_held[[k]]
            k <- k + 1
        }
    }
    amount_held <- bind_rows(list_amount_held)
    amount_held$problem <- factor(amount_held$problem, levels=names(problems))
    amount_held$solution <- factor(amount_held$solution, levels=names(results))
    
    ggplot(amount_held) +
        geom_violin(aes(x=problem, y=relative_held), fill="black") +
        theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
        geom_hline(yintercept=0.15,linetype="dotted") +
        ggtitle(paste(species,"- solutions of",names(problems)[[i.prob]])) ->
        plots[[i.prob]]
}
save(plots,file=paste0(getwd(),"/Results/Plots_amount_held_",species,".RData"))
load(paste0(getwd(),"/Results/Plots_amount_held_",species,".RData"))
png(paste0("Amount_held_prioritizr_",species,".png"),width=15,height=10,units="cm",res=300)
plots[[2]]
dev.off()

# Supplementary Figure: all solutions
lay <- rbind(c(1,2,3),
             c(4,5,6),
             c(7,8,9),
             c(10,11,12),
             c(10,11,12))
for (i in 1 : 9) plots[[i]] <- plots[[i]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
png(paste0("Amount_held_prioritizr_",species,"_SUPPLEMENTARY.png"),width=35,height=30,units="cm",res=300)
marrangeGrob(plots,layout_matrix = lay)
dev.off()


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
            eval_n_summary(select(results[[i.prob]],paste0("solution_",i.sol))) %>% 
            pull(n) ->
            cost[i.prob,i.sol]
    }
}
save(cost,file=paste0(getwd(),"/Results/Cost_prioritizr_",species,".RData"))
load(paste0(getwd(),"/Results/Cost_prioritizr_",species,".RData"))
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


## Define method and num_classes column
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
