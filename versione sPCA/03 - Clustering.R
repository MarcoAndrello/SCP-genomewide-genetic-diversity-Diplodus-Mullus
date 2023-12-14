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
num_axes <- 14
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
save(hopkins, clusternum, file="Results/Diplodus_clustering.RData")
# PAM clustering
pam.res[[1]] <- pam(species_coord, 3,  metric = "euclidean", stand = FALSE)
pam.res[[2]] <- pam(species_coord, 26,  metric = "euclidean", stand = FALSE)
pam.res[[3]] <- pam(species_coord, 2,  metric = "euclidean", stand = FALSE)
# Silhouette plot
png("Figures/Diplodus_silhouette_k3.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[1]], palette = "jco", ggtheme = theme_classic())
dev.off()
png("Figures/Diplodus_silhouette_k26.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[2]], ggtheme = theme_classic())
dev.off()
png("Figures/Diplodus_silhouette_k2.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[3]], palette = "jco", ggtheme = theme_classic())
dev.off()
# Map
pus %>% filter(Diplodus_sargus == 1) %>% mutate(cl = factor(pam.res[[3]]$clustering)) %>%
    select(cl) %>% plot(border=NA)
## write code to make maps of clusters, and show results of the find number of clusters as a table



# Mullus surmuletus
num_axes <- 21
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
save(hopkins, clusternum, file="Results/Mullus_clustering.RData")
# PAM clustering
pam.res <- list()
pam.res[[1]] <- pam(species_coord, 3,  metric = "euclidean", stand = FALSE)
pam.res[[2]] <- pam(species_coord, 30,  metric = "euclidean", stand = FALSE)
pam.res[[3]] <- pam(species_coord, 8,  metric = "euclidean", stand = FALSE)
# Silhouette plot
png("Figures/Mullus_silhouette_k3.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[1]], palette = "jco", ggtheme = theme_classic())
dev.off()
png("Figures/Mullus_silhouette_k30.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[2]], ggtheme = theme_classic())
dev.off()
png("Figures/Mullus_silhouette_k8.png",width=20,height=10,units="cm",res=300)
fviz_silhouette(pam.res[[3]], palette = "jco", ggtheme = theme_classic())
dev.off()
# Map
pus %>% filter(Mullus_surmuletus == 1) %>% mutate(cl = factor(pam.res[[3]]$clustering)) %>%
    select(cl) %>% plot(border=NA)
## write code to make maps of clusters, and show results of the find number of clusters as a table