# Scenarios raptr
rm(list=ls())

# Here set the species you want to work on: "Diplodus" or "Mullus"
species <- "Mullus"

library(tidyverse)
library(sf)
library(terra)
library(raptr)
library(parallel)

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
## Planning units
load("Planning_units.RData")

## Genetic data
g_rast <- rast(paste0(getwd(),genetic_raster))

# Extract genetic values from raster for the centroid of each PU
g_rast_values <- terra::extract(g_rast,pus_centroid)

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
                      proportion = c(0.15,0.95)
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

# Define "gold standard" problem 
ro <- RapUnreliableOpts(BLM=0)
prob_gs <- RapUnsolved(ro, rap_data)
maximum.targets(prob_gs)

# Evaluate the current MR network
selections <- which(rap_data@pu$status==2)
prob_gs_current_MR <- update(prob_gs, b = selections)
space.held(prob_gs_current_MR , y=NULL)
amount.held(prob_gs_current_MR , y=NULL)
spp.plot(prob_gs_current_MR , species=1, y=2)
pus %>% st_drop_geometry %>% filter(status == 2) %>% pull(cost) %>% sum

# Solve the gold standard problem
res_gs <- solve(prob_gs, Threads = threads, verbose=T, NumericFocus= 3L, MIPGap=0.02, NumberSolutions=1L)

save(prob_gs,
     res_gs,
     file=paste0(getwd(),"/Results/Results_raptr_",species,"_95.RData"))

space.held(res_gs, y=NULL)
amount.held(res_gs, y=NULL)
spp.plot(res_gs, species=1, y=2) #, pu.color.palette = c("grey30", "green", "red", "red"))
