# 01 Perfrom spatial PCA and interpolates the results

rm(list=ls())

library(tidyverse)
library(raster)
library(sf)
library(vcfR)
library(adegenet)
library(gstat)

current.dir <- getwd()

# Load planning units
load("Planning_units.RData")

# Create raster covering the extent of pus (for interpolation)
#pusr <- raster(pus,resolution=c(25000,25000))
pusr1 <- raster::raster(pus,resolution=c(10000,10000),crs=st_crs(pus))
pusr2 <- terra::rast(terra::vect(pus), resolution=c(10000,10000))

# Load sampling data: 
# for each cell
# for each Diplodus individuals
# for each Mullus individual
load(paste0(getwd(),"/data/sampling.RData"))


# Diplodus sargus
  
# NEUTRAL LOCI
setwd(paste0(current.dir,"/data"))
n_vcf <- read.vcfR("dip_neutral_7655.vcf.gz")
setwd(current.dir)

# Performing a spatial Pca on the full set of loci is too heavy
# follow the procedure by Jombart
# https://github.com/thibautjombart/adegenet/issues/95

# First do a glPca to retain the most contributing loci
x <- vcfR2genlight(n_vcf)
pca <- glPca(x, nf=2)
toKeep <- loadingplot(pca$loadings^2)$var.idx
save(toKeep, file="Diplodus_toKeep.RData")

# Then perform spatial PCA on the reduced dataset
# read vcf as genind
x <- vcfR2genind(n_vcf)
# Add population information
x$pop <- Diplodus_sampling$SamplingCell
# Convert to genpop
xpop <- genind2genpop(x)
# Match rownames of xpop to sampling cells in cell_sampling to retrieve geographic coordinates
match_xpop_to_cellsampling <- match(row.names(xpop$tab),cell_sampling$SamplingCell)
# Retrieve geographic coordinates
Diplodus_coord <- cell_sampling[match_xpop_to_cellsampling,c("Longitude","Latitude")]
rm(match_xpop_to_cellsampling)
save(Diplodus_coord,file="Diplodus_coord.RData")
# Add spatial coordinates
other(xpop)$xy <- Diplodus_coord 
names(other(xpop)$xy) <- c("x","y")
# Keep only the most contributing loci
load("Diplodus_toKeep.RData")
xpop <- xpop[,loc=toKeep]
# Perform spatial PCA
spca_neutr <- spca(xpop)
save(spca_neutr,file=paste0(getwd(),"/Results sPCA/Diplodus_spca_neutr.RData"))

# Interpolation
# might want to read this:
# http://132.72.155.230:3838/r/spatial-interpolation-of-point-data.html
# and see why terra::interpolate does not work with gstat on sf objects
# https://github.com/rspatial/terra/issues/208
# (we need terra::interpolate and not raster::interpolate because only
# terra::rast seems to deal well with our CRS (Albers)
# trying a spatVect
load(paste0(getwd(),"/Results sPCA/Diplodus_spca_neutr.RData"))
obs.data <- data.frame(spca1=spca_neutr$ls[,1],
                       spca2=spca_neutr$ls[,2],
                       lon=Diplodus_coord$Longitude,
                       lat=Diplodus_coord$Latitude)
obs <- vect(obs.data,
                geom=c("lon","lat"),
                crs="EPSG:4326")
obs <- project(obs,"ESRI:102013")
gs1 <- gstat(formula=spca1~1, locations=obs, nmax=5, set=list(idp=0))
gs1 <- gstat(formula=spca1~1, obs, locations=~x+y, nmax=5, set=list(idp=0))
ng1 <- terra::interpolate(pusr2, gs1)
gOK2 <- gstat(NULL, "log.zinc", log(zinc)~1, locations = meuse_sf, model = mv2)
gOK <- gstat(NULL, "log.zinc", log(zinc)~1, meuse, locations=~x+y, model=mv)





load(paste0(getwd(),"/Results sPCA/Diplodus_spca_neutr.RData"))
# Create an sf object with observed values
obs <- st_as_sf(data.frame(spca1=spca_neutr$ls[,1],
                           spca2=spca_neutr$ls[,2],
                           lon=Diplodus_coord$Longitude,
                           lat=Diplodus_coord$Latitude),
         coords=c("lon","lat"),
         crs=4326)
# Transform to Albers equal area conic
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
# Define nearest neighbor model and interpolate on pusr

## with a model built with an `sf` object you need to provide custom function
gs1 <- gstat(formula=spca~1, data=obs, nmax=5, set=list(idp = 0))

interpolate_gstat <- function(model, x, crs, ...) {
    v <- st_as_sf(x, coords=c("x", "y"), crs=crs)
    p <- predict(model, v, ...)
    as.data.frame(p)[,1:2]
}

ng1 <- interpolate(pusr2, gs1, debug.level=0, fun=interpolate_gstat, crs=crs(pusr2), index=1)

# NON FUNZIONA.
# E COMUNQUE, RIJ_MATRIX NON FUNZIONEREBBE CON SPATRASTER. BISOGNA FARLO FUNZIONARE SUL RASTER E CON LA CRS GIUSTA

###


## First neutral sPCA axis
gs1 <- gstat(formula=spca1~1, locations=obs, nmax=5, set=list(idp=0))
# ng1 <- interpolate(pusr, gs1)
ng1 <- interpolate(pusr1, gs1)
ng1 <- terra::interpolate(rast(pusr1), gs1)
## Second neutral sPCA axis
gs2 <- gstat(formula=spca2~1, locations=obs, nmax=5, set=list(idp=0))
ng2 <- interpolate(pusr, gs2)
plot(ng1)
plot(ng2)
Diplodus_ng <- brick(ng1, ng2)
names(Diplodus_ng) <- c("Neutral_sPCA_Axis_1", "Neutral_sPCA_Axis_2")
save(Diplodus_ng,file=paste0(getwd,"/Results sPCA/Diplodus_neutral_SPCA.RData"))
#writeRaster(Diplodus_ng,filename="Diplodus_neutral_SPCA.tif")
rm(obs,gs1,gs2,ng1,ng2,Diplodus_ng)
# https://rspatial.org/raster/analysis/4-interpolation.html


#### OUTLIER LOCI
rm(n_vcf, x, xpop, toKeep)
setwd(paste0(current.dir,"/data"))
n_vcf <- read.vcfR("dip_adaptive_413.vcf.gz")
setwd(current.dir)
x <- vcfR2genind(n_vcf)
x$pop <- Diplodus_sampling$SamplingCell
xpop <- genind2genpop(x)
other(xpop)$xy <- Diplodus_coord
names(other(xpop)$xy) <- c("x","y")
spca_outlier <- spca(xpop)
save(spca_outlier,file=paste0(getwd(),"/Results sPCA/Diplodus_spca_outlier.RData"))

# Nearest neighbor interpolation
load(paste0(getwd(),"/Results sPCA/Diplodus_spca_outlier.RData"))
# Create an sf object with observed values
obs <- st_as_sf(data.frame(spca1=spca_outlier$ls[,1],
                           spca2=spca_outlier$ls[,2],
                           lon=Diplodus_coord$Longitude,
                           lat=Diplodus_coord$Latitude),
                coords=c("lon","lat"),
                crs=4326)
# Transform to Albers equal area conic
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
gs1 <- gstat(formula=spca1~1, locations=obs, nmax=5, set=list(idp=0))
og1 <- interpolate(pusr,gs1)
gs2 <- gstat(formula=spca2~1, locations=obs, nmax=5, set=list(idp=0))
og2 <- interpolate(pusr,gs2)
Diplodus_og <- brick(og1, og2)
names(Diplodus_og) <- c("Outlier_sPCA_Axis_1", "Outlier_sPCA_Axis_2")
save(Diplodus_og,file=paste0(getwd(),"/Results sPCA/Diplodus_outlier_SPCA.RData"))
# writeRaster(Diplodus_og,filename="Diplodus_outlier_SPCA.tif")
rm(obs,gs1,gs2,og1,og2,Diplodus_og)
rm(n_vcf, x, xpop)
rm(spca_neutr, spca_outlier, Diplodus_coord)


###

# Mullus surmuletus

# NEUTRAL LOCI
# Follows the exact same code for Diplodus sargus above, hence it is not commented

setwd(paste0(current.dir,"/data"))
n_vcf <- read.vcfR("mul_neutral_2462.vcf.gz")
setwd(current.dir)

x <- vcfR2genlight(n_vcf)
pca <- glPca(x, nf=2)
toKeep <- loadingplot(pca$loadings^2)$var.idx
save(toKeep, file="Mullus_toKeep.RData")

x <- vcfR2genind(n_vcf)
x$pop <- Mullus_sampling$SamplingCell
xpop <- genind2genpop(x)
match_xpop_to_cellsampling <- match(row.names(xpop$tab),cell_sampling$SamplingCell)
Mullus_coord <- cell_sampling[match_xpop_to_cellsampling,c("Longitude","Latitude")]
rm(match_xpop_to_cellsampling)
other(xpop)$xy <- Mullus_coord 
names(other(xpop)$xy) <- c("x","y")
load("Mullus_toKeep.RData")
xpop <- xpop[,loc=toKeep]
spca_neutr <- spca(xpop)
save(spca_neutr,file=paste0(getwd(),"/Results sPCA/Mullus_spca_neutr.RData"))

load(paste0(getwd(),"/Results sPCA/Mullus_spca_neutr.RData"))
obs <- st_as_sf(data.frame(spca1=spca_neutr$ls[,1],
                           spca2=spca_neutr$ls[,2],
                           lon=Mullus_coord$Longitude,
                           lat=Mullus_coord$Latitude),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
gs1 <- gstat(formula=spca1~1, locations=obs, nmax=5, set=list(idp=0))
ng1 <- interpolate(pusr, gs1)
gs2 <- gstat(formula=spca2~1, locations=obs, nmax=5, set=list(idp=0))
ng2 <- interpolate(pusr, gs2)
plot(ng1)
plot(ng2)
Mullus_ng <- brick(ng1, ng2)
names(Mullus_ng) <- c("Neutral_sPCA_Axis_1", "Neutral_sPCA_Axis_2")
save(Mullus_ng,file=paste0(getwd(),"/Results sPCA/Mullus_neutral_SPCA.RData"))
# writeRaster(Mullus_ng,filename="Mullus_neutral_SPCA.tif") ### DOES NOT WORK
rm(obs,gs1,gs2,ng1,ng2,Mullus_ng)

#### OUTLIER LOCI
rm(n_vcf, x, xpop, toKeep)
setwd(paste0(current.dir,"/data"))
n_vcf <- read.vcfR("mul_adaptive_291.vcf.gz")
setwd(current.dir)
x <- vcfR2genind(n_vcf)
x$pop <- Mullus_sampling$SamplingCell
xpop <- genind2genpop(x)
other(xpop)$xy <- Mullus_coord
names(other(xpop)$xy) <- c("x","y")
spca_outlier <- spca(xpop)
save(spca_outlier,file=paste0(getwd(),"/Results sPCA/Mullus_spca_outlier.RData"))

load(paste0(getwd(),"/Results sPCA/Mullus_spca_outlier.RData"))
obs <- st_as_sf(data.frame(spca1=spca_outlier$ls[,1],
                           spca2=spca_outlier$ls[,2],
                           lon=Mullus_coord$Longitude,
                           lat=Mullus_coord$Latitude),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
gs1 <- gstat(formula=spca1~1, locations=obs, nmax=5, set=list(idp=0))
og1 <- interpolate(pusr,gs1)
gs2 <- gstat(formula=spca2~1, locations=obs, nmax=5, set=list(idp=0))
og2 <- interpolate(pusr,gs2)
Mullus_og <- brick(og1, og2)
names(Mullus_og) <- c("Outlier_sPCA_Axis_1", "Outlier_sPCA_Axis_2")
save(Mullus_og,file=paste0(getwd(),"/Results sPCA/Mullus_outlier_SPCA.RData"))
# writeRaster(Mullus_og,filename="Mullus_outlier_SPCA.tif") # does not work
rm(obs,gs1,gs2,og1,og2,Mullus_og)
rm(n_vcf, x, xpop)

rm(spca_neutr, spca_outlier, Mullus_coord)

# FIGURES
load(paste0(getwd(),"/Results sPCA/Diplodus_neutral_SPCA.RData"))
load(paste0(getwd(),"/Results sPCA/Diplodus_outlier_SPCA.RData"))
load(paste0(getwd(),"/Results sPCA/Mullus_neutral_SPCA.RData"))
load(paste0(getwd(),"/Results sPCA/Mullus_outlier_SPCA.RData"))

st_read(paste0(getwd(),"/../../Maps/ne_110m_land.shp")) %>%
st_transform(st_crs(Diplodus_ng[[1]])) %>%
  st_geometry() %>% 
  st_crop(st_bbox(Diplodus_ng[[1]])) ->
  countries

png("Diplodus_ng_Axis1.png",width=16,height=9.4,units="cm", res=300)
plot(Diplodus_ng[[1]],main="D. sargus, neutral, Axis 1")
plot(countries,add=T)
dev.off()

png("Diplodus_ng_Axis2.png",width=16,height=9.4,units="cm", res=300)
plot(Diplodus_ng[[2]],main="D. sargus, neutral, Axis 2")
plot(countries,add=T)
dev.off()

png("Diplodus_og_Axis1.png",width=16,height=9.4,units="cm", res=300)
plot(Diplodus_og[[1]],main="D. sargus, outlier, Axis 1")
plot(countries,add=T)
dev.off()

png("Diplodus_og_Axis2.png",width=16,height=9.4,units="cm", res=300)
plot(Diplodus_og[[2]],main="D. sargus, outlier, Axis 2")
plot(countries,add=T)
dev.off()


png("Mullus_ng_Axis1.png",width=16,height=9.4,units="cm", res=300)
plot(Mullus_ng[[1]],main="M. surmuletus, neutral, Axis 1")
plot(countries,add=T)
dev.off()

png("Mullus_ng_Axis2.png",width=16,height=9.4,units="cm", res=300)
plot(Mullus_ng[[2]],main="M. surmuletus, neutral, Axis 2")
plot(countries,add=T)
dev.off()

png("Mullus_og_Axis1.png",width=16,height=9.4,units="cm", res=300)
plot(Mullus_og[[1]],main="M. surmuletus, outlier, Axis 1")
plot(countries,add=T)
dev.off()

png("Mullus_og_Axis2.png",width=16,height=9.4,units="cm", res=300)
plot(Mullus_og[[2]],main="M. surmuletus, outlier, Axis 2")
plot(countries,add=T)
dev.off()

# # Some more plots
# plot(spca_neutr$ls[,1], spca_neutr$ls[,2], xlab="Neutral loci sPCA axis 1", ylab="Neutral loci sPCA axis 2")
# plot(spca_outlier$ls[,1], spca_outlier$ls[,2], xlab="Adaptive loci sPCA axis 1", ylab="Adaptive loci sPCA axis 2")
# plot(spca_neutr$ls[,1], spca_outlier$ls[,1], xlab="Neutral loci sPCA axis 1", ylab="Adaptive loci sPCA axis 1")
# 
# # Plot graphs with longitude
# plot(Diplodus_coord$Longitude, spca_neutr$ls[,1], xlab="Longitude", ylab="First sPCA axis (neutral loci)")
# plot(Diplodus_coord$Longitude, spca_outlier$ls[,1], xlab="Longitude", ylab="First sPCA axis (adaptive loci)")
