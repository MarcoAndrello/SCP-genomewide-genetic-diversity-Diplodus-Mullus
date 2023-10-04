# 02 Perform PCA and interpolates the results

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(adegenet)
library(gstat)
library(rnaturalearth)

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

# Load merged (neutral + adaptive) dataset (but this exlcudes the non-outlier that are  not in HWE)
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

# Perform (spatial) PCA
# spca_neutr <- spca(xpop)
# save(spca_neutr,file=paste0(getwd(),"/Results sPCA/Diplodus_spca_neutr.RData"))
# NEED TO DO IT ONLY ON LOCI, NOT ALLELES
pca <- dudi.pca(xpop)
save(pca,file=paste0(getwd(),"/Results May_2023/Diplodus_pca_pop.RData"))

# Plot PCA (biplot and screeplot)
load(paste0(getwd(),"/Results May_2023/Diplodus_pca_pop.RData"))
png(paste0("PCA_biplot_Diplodus.png"),width=15,height=15,units="cm",res=300)
plot(pca$li[,1],pca$li[,2],pch=16,xlab="Axis 1",ylab="Axis 2",col="gray", main = "PCA Diplodus sargus")
dev.off()
screeplot(pca)
pca$eig/sum(pca$eig)

 # Interpolation
# load(paste0(getwd(),"/Results May_2023/Diplodus_pca.RData"))
# Define obs sf object so we can convert the lat-lon coordinates to Albers equal area
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
Diplodus_coord_Albers <- st_coordinates(obs)
rm(obs)
# Loop on PCA axes
pca_interp <- list()
for (i.axis in 1 :ncol(pca$li)){
    if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=pca$li[,i.axis], X=Diplodus_coord_Albers[,1], Y=Diplodus_coord_Albers[,2])
    # gstat model
    gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0))
    # Interpolation
    a <- terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))
    pca_interp[[i.axis]] <- a[[1]]
}
# Save as raster
Diplodus_pca <- rast(pca_interp)
writeRaster(Diplodus_pca, file="Diplodus_allAxes_8068.grd")

# Plotting the rasters
g_rast <- rast(paste0(getwd(),"/Results May_2023/Diplodus_allAxes_8068.grd"))
names(g_rast) <- paste("Axis",c(1:nlyr(g_rast)))
# Read countries for plotting maps
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(g_rast)) -> countries
png(paste0("Genetic_rasters_Diplodus.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(6,3))
is.last <- F
for (i in 1 : 17) {
    if (i == 17) is.last=T
    plot(g_rast,i,mar=c(0.5,0.1,1,0.1),axes=F,legend=is.last,main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()

# Absolute scaled values averaged over layers
g_rast_scaled_abs <- abs(scale(g_rast))
png(paste0("Average_raster_Diplodus.png"),width=10,height=8,units="cm",res=300)
plot(mean(g_rast_scaled_abs[[1:17]]),main="Diplodus sargus",col=rev(heat.colors(20)),cex.main=0.5)
polys(countries,col="gray",lwd=0.01)
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
# Perform PCA
pca <- dudi.pca(xpop)
save(pca,file=paste0(getwd(),"/Results May_2023/Mullus_pca_pop.RData"))
# Plot PCA (biplot and screeplot)
load(paste0(getwd(),"/Results May_2023/Mullus_pca_pop.RData"))
png(paste0("PCA_biplot_Mullus.png"),width=15,height=15,units="cm",res=300)    
plot(pca$li[,1],pca$li[,2],pch=16,xlab="Axis 1",ylab="Axis 2",col="gray", main = "PCA Mullus surmuletus")
dev.off()
screeplot(pca)
pca$eig/sum(pca$eig)

# Interpolation
# load(paste0(getwd(),"/Results May_2023/Mullus_pca.RData"))
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
Mullus_coord_Albers <- st_coordinates(obs)
rm(obs)
pca_interp <- list()
for (i.axis in 1 : ncol(pca$li)){
  if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
  data <- data.frame(pca=pca$li[,i.axis], X=Mullus_coord_Albers[,1], Y=Mullus_coord_Albers[,2])
  gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0))
  a <- terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))
  pca_interp[[i.axis]] <- a[[1]]
}
Mullus_pca <- rast(pca_interp)
writeRaster(Mullus_pca, file="Mullus_allAxes_2753.grd")

# Plotting the rasters
g_rast <- rast(paste0(getwd(),"/Results May_2023/Mullus_allAxes_2753.grd"))
names(g_rast) <- paste("Axis",c(1:nlyr(g_rast)))
# Read countries for plotting maps
ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(g_rast)) -> countries
is.last <- F
png(paste0("Genetic_rasters_Mullus_1.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(6,3))
for (i in 1 : 18) {
    plot(g_rast,i,mar=c(0.5,0.1,1,0.1),axes=F,legend=F,main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()
png(paste0("Genetic_rasters_Mullus_2.png"),width=20,height=24,units="cm",res=300)
par(mfrow=c(6,3))
for (i in 19 : 26) {
    if (i == 26) is.last=T
    plot(g_rast,i,mar=c(0.5,0.1,1,0.1),axes=F,legend=is.last,main=paste("Axis",i))
    polys(countries,col="gray",lwd=0.01)
}
dev.off()
# Absolute scaled values averaged over layers
g_rast_scaled_abs <- abs(scale(g_rast))
png(paste0("Average_raster_Mullus.png"),width=10,height=8,units="cm",res=300)
plot(mean(g_rast_scaled_abs[[1:26]]),main="Mullus surmuletus",col=rev(heat.colors(20)),cex.main=0.5)
polys(countries,col="gray",lwd=0.01)
dev.off()

xtab <- xpop$tab
xtab <- xtab[,match(levels(xpop@loc.fac),xpop@loc.fac)]
pcal <- dudi.pca(xtab)
plot(pca$li[,3],pcal$li[,4])
plot(pca$li[,4],pcal$li[,3])

# Interpolation code:
  # # Using terra
  # obs1 <- cbind(st_drop_geometry(obs),st_coordinates(obs))
  # gs <- gstat(formula=pca1~1, locations=~X+Y, data=obs1, nmax=5, set=list(idp=0))
  # a1 <-terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))
  # # Using raster
  # gs2 <- gstat(formula=pca1~1, locations=obs, nmax=5, set=list(idp=0))
  # a2 <- raster::interpolate(pusr, gs2)
