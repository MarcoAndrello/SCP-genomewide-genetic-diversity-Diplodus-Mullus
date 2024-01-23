# 09 - Mapping

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(rnaturalearth)
library(tmap)

load("Planning_units.RData")

ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(pus)) -> countries

################################################################################
# Figure S1
# Map of conservation cost
################################################################################
png("Figures/Map_cost.png",
    width=20,height=12,res=600,units="cm")
print(
    tm_shape(mutate(pus, cost=cost/1000)) +
        tm_fill("cost", style="log10_pretty", legend.is.portrait = F, legend.show = T) +
        tm_shape(filter(pus, status_0 == 1))+
        tm_borders(col = "blue", lwd=0.2) +
        tm_shape(countries, bbox = res) +
        tm_polygons(col="lightgray", lwd=0.2) +
        tm_legend(legend.outside = T, legend.outside.position = "bottom") +
        tm_layout(main.title="Conservation cost (M€)", main.title.position = "center")
)
dev.off()


################################################################################
# Figure S2
# Map of sampling points
################################################################################
load("data/sampling.RData")
Diplodus_sampling %>%
    group_by(SamplingCell) %>%
    summarise(Long=mean(Longitude), Lat=mean(Latitude)) %>%
    st_as_sf(coords=c("Long","Lat"),
             crs=st_crs(4326)) %>%
    st_transform(st_crs(pus)) -> Diplodus_points
Mullus_sampling %>%
    group_by(SamplingCell) %>%
    summarise(Long=mean(Longitude), Lat=mean(Latitude)) %>%
    st_as_sf(coords=c("Long","Lat"),
             crs=st_crs(4326)) %>%
    st_transform(st_crs(pus)) -> Mullus_points

png("Figures/Sampling_points_Diplodus.png",width=20,height=8,units="cm",res=600)
tm_shape(mutate(pus, Diplodus_sargus=factor(Diplodus_sargus))) +
    tm_fill("Diplodus_sargus",
            palette = c("grey90","blue"),
            legend.show=F) +
    tm_shape(countries, bbox = pus) +
    tm_polygons(col="lightgray", lwd=0.2) +
    tm_shape(Diplodus_points) + 
    tm_dots(size=0.25)
dev.off()
png("Figures/Sampling_points_Mullus.png",width=20,height=8,units="cm",res=600)
tm_shape(mutate(pus, Mullus_surmuletus=factor(Mullus_surmuletus))) +
    tm_fill("Mullus_surmuletus",
            palette = c("grey90","blue"),
            legend.show=F) +
    tm_shape(countries, bbox = pus) +
    tm_polygons(col="lightgray", lwd=0.2) +
    tm_shape(Mullus_points) + 
    tm_dots(size=0.25)
dev.off()


################################################################################
# Figure S6 and S7
# Spatially interpolated PCA scores for D. sargus and M. surmuletus
################################################################################

# Diplodus
g_rast <- rast("Results/Diplodus_allAxes_8068.grd")
png(paste0("Genetic_rasters_Diplodus_1.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(4,3))
for (i in 1 : 12) {
    plot(g_rast,i,mar=c(2,0.5,0.5,0.5),axes=F,legend="bottom",main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()
png(paste0("Genetic_rasters_Diplodus_2.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(4,3))
for (i in 13 : 17) {
    plot(g_rast,i,mar=c(2,0.5,0.5,0.5),axes=F,legend="bottom",main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()

# Mullus
g_rast <- rast("Results/Mullus_allAxes_2753.grd")
# Read countries for plotting maps
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(g_rast)) -> countries
png(paste0("Genetic_rasters_Mullus_1.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(4,3))
for (i in 1 : 12) {
    plot(g_rast,i,mar=c(2,0.5,0.5,0.5),axes=F,legend="bottom",main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()
png(paste0("Genetic_rasters_Mullus_2.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(4,3))
for (i in 13 : 24) {
    plot(g_rast,i,mar=c(2,0.5,0.5,0.5),axes=F,legend="bottom",main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()
png(paste0("Genetic_rasters_Mullus_3.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(4,3))
for (i in 25 : 26) {
    plot(g_rast,i,mar=c(2,0.5,0.5,0.5),axes=F,legend="bottom",main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()


################################################################################
# Figure S8
# Map of priority sites for protecting 15% of species ranges without explicit genetic objectives. 
################################################################################
load("Results/Results_Extension_Amount.RData")

# Calculate selection frequency
selections <- res_so %>%
    st_drop_geometry() %>% select(paste0("solution_",c(1:100)))
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

png("Figures/Maps_selections/Map_selection_frequency_speciesOnly.png",
    width=20,height=12,res=600,units="cm")
print(
    tm_shape(pus_map) +
        tm_fill("selection_frequency", palette=c("gold",brewer.pal(5,"Greens")),
                legend.is.portrait = F, legend.show = T) +
        tm_shape(countries, bbox = res) +
        tm_polygons(col="lightgray", lwd=0.2) +
        tm_legend(legend.outside = T, legend.outside.position = "bottom") +
        tm_layout(main.title="Range only", main.title.position = "center")
)
dev.off()
