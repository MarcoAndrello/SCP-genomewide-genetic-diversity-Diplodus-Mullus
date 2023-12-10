# phylin
library(phylin)
library(adegenet)
library(tidyverse)
library(sf)
load("Data_for_PCA_Diplodus.RData")
d.gen <- dist.genpop(xpop)
xpop@other$xy %>% st_as_sf(coords=c("x","y"),crs=4326) %>% st_distance -> r.dist
units::drop_units(r.dist)/1000 -> r.dist
gv <- gen.variogram(r.dist, d.gen)
plot(gv)
gv$n
# but then, to interpolate, we need presence or absence of lineages. we do not have such a thing with SNPs (all are present everywhere)