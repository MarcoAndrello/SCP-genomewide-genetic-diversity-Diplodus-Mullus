# spBayes
library(spBayes)
library(adegenet)
library(tidyverse)
library(sf)
library(terra)
load("Data_for_PCA_Diplodus.RData")
load("Planning_units.RData")
pus %>% st_transform(st_crs(4326)) -> pus
xpop@tab[,1] -> y
xpop@tab[,1] + xpop@tab[,2] -> N
xpop@other$xy -> coords
knots <- expand.grid(x=seq(-6,36,0.5),y=seq(30,45,0.5))
# pust <- rast(pus, resolution = c(0.5, 0.5), crs = "EPSG:4326")
m <- spGLM(y~1, family="binomial", weights=N, coords=as.matrix(coords), knots=as.matrix(knots),
           n.samples=1000, cov.model="exponential")

data <- st_as_sf(data.frame(counts=y/N,coords),coords=c("x","y"),crs=st_crs(4326))
plot(data,pch=16)
