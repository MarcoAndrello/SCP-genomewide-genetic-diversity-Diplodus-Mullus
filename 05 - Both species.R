# For both species, evaluate the space held by the current MR network

rm(list=ls())

species <- "Mullus"

library(tidyverse)
library(sf)
library(terra)
library(prioritizr)
library(raptr)

if (species == "Diplodus") {
    name_species = "Diplodus_sargus"
    genetic_raster = "/Results/Diplodus_allAxes_8068.grd"
    name_species_feat = "Diplodus_a"
    num_axes = 17
}
if (species == "Mullus") {
    name_species = "Mullus_surmuletus"
    genetic_raster = "/Results/Mullus_allAxes_2753.grd"
    name_species_feat = "Mullus_a"
    num_axes = 26
}


# Number of assigned threads
threads <- 8L #as.integer(max(1, detectCores() - 2))

# Load data
# Planning units
load("Planning_units.RData")

## Genetic data
g_rast <- rast(paste0(getwd(),genetic_raster))
g_rast_values <- terra::extract(g_rast,pus_centroid)

## mi fermo qui

# Create an object containing the coordinates of the occupied PU
# in all genetic spaces (i.e. for all genetic PCA axes)(need this for planning.unit.points)
pus_g_values <- cbind(g_rast_values, st_drop_geometry(pus))
names(pus_g_values)[1:ncol(g_rast_values)] <- paste0("pca",sprintf("%02d",1:ncol(g_rast_values)))
species_coord <- filter(pus_g_values, !!as.name(name_species) == 1) 
species_coord <- species_coord[1:num_axes]
rm(pus_g_values)

# Create species probabilities
pu.species.probabilities <- data.frame(species=1L,
                                       pu=c(1:nrow(pus)),
                                       value=pus[[name_species]])
pu.species.probabilities <- filter(pu.species.probabilities, value==1)

# Convert to PolySet and calculate boundary
pus %>%  convert2PolySet() -> polygons
boundary <- raptr::calcBoundaryData(polygons)

# Many PUs have the same multidimensional coordinates in the PCA space (this is due to the interpolation method)
# Count number of PUs per multidimensional coordinate
species_coord %>% group_by_all() %>% count -> coords
# The number of PUs per coordinate will be the weight of that coordinate
coords[,num_axes+1] %>% pull(n) -> weights
weights <- weights/sum(weights)
# Remove last column (n, the counts)
coords <- coords[,1:num_axes]
# Make demand points
dp <- DemandPoints(as.matrix(coords), weights)

# Attribute space       # In this list there will be AttributeSpace for each species
genetic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=as.matrix(species_coord),
                                                                                ids=pu.species.probabilities$pu), # only occupied PUs
                                      demand.points = dp,
                                      species=1L))
# Attribute spaces for this axis (both species)
attribute.spaces <- list()
attribute.spaces[[1]] <-
    AttributeSpaces(
        spaces = genetic_spaces,
        name = "genetic"
    )
# Create targets
## Geographic targets to 0
targets <- data.frame(species = 1L,
                      target = as.integer(0:1),
                      proportion = c(0.15,0.75)
)

pus$status[pus$status==1] <- pus$status[pus$status==1]+1
pus$status <- as.integer(pus$status)

rap_data <- RapData(pu = st_drop_geometry(pus),
                    species = data.frame(name=species),
                    targets = targets,
                    pu.species.probabilities = pu.species.probabilities,
                    attribute.spaces = attribute.spaces,
                    boundary = boundary,
                    polygons = polygons)

# Define "gold standard" problem and solve it
ro <- RapUnreliableOpts(BLM=0)
prob_gs <- RapUnsolved(ro, rap_data)
maximum.targets(prob_gs)
res_gs <- solve(prob_gs, Threads = threads, verbose=T, NumericFocus= 3L, MIPGap=0.02, NumberSolutions=100L)




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

