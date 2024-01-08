# 09 - Mapping

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(rnaturalearth)
library(tmap)

load("Planning_units.RData")

ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(pus)) -> countries

# Map of sampling points
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