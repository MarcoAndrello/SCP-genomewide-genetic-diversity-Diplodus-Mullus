# Clustering

rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
library(cluster)
library(NbClust)
library(factoextra)
library(hopkins)

# Load data
## Planning units
load("Planning_units.RData")

# Open genetic raster and Extract values for the centroid of each PU
g_rast_values <- list()
g_rast <- rast(paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))
g_rast_values[[1]] <- terra::extract(g_rast,pus_centroid)
g_rast <- rast(paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))
g_rast_values[[2]] <- terra::extract(g_rast,pus_centroid)
rm(g_rast)


# https://rstudio-pubs-static.s3.amazonaws.com/375287_5021917f670c435bb0458af333716136.html
# Diplodus sargus
num_axes <- 17
# Create an object containing the coordinates of the occupied PU
# in all genetic PCA axes: need this for finding clusters
pus_g_values <- cbind(g_rast_values[[1]], st_drop_geometry(pus))
names(pus_g_values)[1:ncol(g_rast_values[[1]])] <- paste0("spca",sprintf("%02d",1:ncol(g_rast_values[[1]])))
species_coord <- filter(pus_g_values, Diplodus_sargus == 1)
species_coord <- species_coord[1:num_axes]
rm(pus_g_values)
# Hopkins test to see if the data can be clustered
hopkins <- hopkins::hopkins(species_coord,m=500)
# pca to get a flavour of the possible number of clusters
fviz_pca_ind(prcomp(species_coord), title = "PCA - diplodus", palette = "jco", geom = "point", ggtheme = theme_classic(), legend = "bottom")
# run different algorithms to find the number of clusters
clusternum <- NbClust(species_coord, distance="euclidean", method="kmeans", index="all", max.nc = 30)
# print results
clusternum$Best.nc[1,] %>% table %>% sort(decreasing=T) %>% names %>% as.numeric %>% print
save(hopkins, clusternum, file="Results/Clustering_Diplodus.RData")


# Mullus surmuletus
num_axes <- 26
# Create an object containing the coordinates of the occupied PU
# in all genetic PCA axes: need this for finding clusters
pus_g_values <- cbind(g_rast_values[[2]], st_drop_geometry(pus))
names(pus_g_values)[1:ncol(g_rast_values[[2]])] <- paste0("spca",sprintf("%02d",1:ncol(g_rast_values[[2]])))
species_coord <- filter(pus_g_values, Mullus_surmuletus == 1)
species_coord <- species_coord[1:num_axes]
rm(pus_g_values)
# Hopkins test to see if the data can be clustered
hopkins <- hopkins::hopkins(species_coord,m=500)
# pca to get a flavour of the possible number of clusters
fviz_pca_ind(prcomp(species_coord), title = "PCA - Mullus", palette = "jco", geom = "point", ggtheme = theme_classic(), legend = "bottom")
# run different algorithms to find the number of clusters
clusternum <- NbClust(species_coord, distance="euclidean", method="kmeans", index="all", max.nc = 30)
# print results
clusternum$Best.nc[1,] %>% table %>% sort(decreasing=T) %>% names %>% as.numeric %>% print
save(hopkins, clusternum, file="Results/Clustering_Mullus.RData")