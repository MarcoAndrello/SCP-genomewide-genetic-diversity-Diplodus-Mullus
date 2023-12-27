# 02 Perform PCA and interpolates the results

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(adegenet)
library(gstat)
library(rnaturalearth)
library(tictoc)


current.dir <- getwd()

# Load planning units
load("Planning_units.RData")

# Create raster covering the extent of pus (for interpolation)
pust <- rast(pus, resolution = c(10000, 10000), crs = crs(pus))

# Load sampling data: 
# for each cell
# for each Diplodus individuals
# for each Mullus individual
load(paste0(getwd(),"/data/sampling.RData"))



####################
# Diplodus sargus
####################

# Load merged (neutral + adaptive) dataset (but this excludes the non-outlier that are not in HWE)
load(paste0(getwd(),"/data/Diplodus_8068.RData"))
x <- Diplodus_8068; rm(Diplodus_8068)
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
save(xpop,file="Data_for_PCA_Diplodus.RData")

# Perform PCA
load("Data_for_PCA_Diplodus.RData")
pca <- dudi.pca(xpop)
save(pca,file=paste0(getwd(),"/Results/Diplodus_pca_pop.RData"))

# Plot PCA (biplot and screeplot)
load(paste0(getwd(),"/Results/Diplodus_pca_pop.RData"))
png(paste0("PCA_biplot_Diplodus.png"),width=15,height=15,units="cm",res=300)
plot(pca$li[,1],pca$li[,2],pch=16,xlab="Axis 1",ylab="Axis 2",col="gray", main = "PCA Diplodus sargus")
dev.off()
screeplot(pca)
pca$eig/sum(pca$eig)

# Interpolation
load("Data_for_PCA_Diplodus.RData")
load("Results/Diplodus_pca_pop.RData")
# Define obs sf object so we can convert the lat-lon coordinates to Albers equal area
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
coord_Albers <- st_coordinates(obs)
rm(obs)
# Loop on PCA axes
pca_interp <- list()
for (i.axis in 1 : 17){
    if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=pca$li[,i.axis], X=coord_Albers[,1], Y=coord_Albers[,2])
    # gstat model
    # gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0)) ## nn
    gs <- gstat(formula=pca~1, locations=~X+Y, data=data) ## idw
    # Interpolation
    terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))$var1.pred %>%
        mask(filter(pus, Diplodus_sargus==1)) ->
    pca_interp[[i.axis]]
}

# Save as raster
names(pca_interp) <- paste("Axis",c(1:17))
pca_rast <- rast(pca_interp)
writeRaster(pca_rast, file=paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))

# Plotting the rasters - Figure S3
g_rast <- rast(paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))
# Read countries for plotting maps
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(g_rast)) -> countries
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



####################
# Mullus surmuletus
####################
load(paste0(getwd(),"/data/Mullus_2753.RData"))
x <- Mullus_2753; rm(Mullus_2753)
x$pop <- Mullus_sampling$SamplingCell
xpop <- genind2genpop(x)
match_xpop_to_cellsampling <- match(row.names(xpop$tab),cell_sampling$SamplingCell)
Mullus_coord <- cell_sampling[match_xpop_to_cellsampling,c("Longitude","Latitude")]
rm(match_xpop_to_cellsampling)
save(Mullus_coord,file="Mullus_coord.RData")
other(xpop)$xy <- Mullus_coord 
names(other(xpop)$xy) <- c("x","y")
save(xpop,file="Data_for_PCA_Mullus.RData")

# Perform PCA
pca <- dudi.pca(xpop)
save(pca,file=paste0(getwd(),"/Results/Mullus_pca_pop.RData"))
# Plot PCA (biplot and screeplot)
load(paste0(getwd(),"/Results/Mullus_pca_pop.RData"))
png(paste0("PCA_biplot_Mullus.png"),width=15,height=15,units="cm",res=300)    
plot(pca$li[,1],pca$li[,2],pch=16,xlab="Axis 1",ylab="Axis 2",col="gray", main = "PCA Mullus surmuletus")
dev.off()
screeplot(pca)
pca$eig/sum(pca$eig)

# Interpolation
load("Data_for_PCA_Mullus.RData")
load("Results/Mullus_pca_pop.RData")
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
coord_Albers <- st_coordinates(obs)
rm(obs)
pca_interp <- list()
for (i.axis in 1 : 26){
    if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=pca$li[,i.axis], X=coord_Albers[,1], Y=coord_Albers[,2])
    # gstat model
    # gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0)) ## nn
    gs <- gstat(formula=pca~1, locations=~X+Y, data=data) ## idw
    # Interpolation
    terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))$var1.pred %>%
        mask(filter(pus, Mullus_surmuletus==1)) ->
        pca_interp[[i.axis]]
}

# Save as raster
names(pca_interp) <- paste("Axis",c(1:26))
pca_rast <- rast(pca_interp)
writeRaster(pca_rast, file=paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))


# Plotting the rasters - Figure S4
g_rast <- rast(paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))
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

