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

# Try sPCA
load("Data_for_PCA_Diplodus.RData")
xpop1 <- xpop[,loc=c(1:100)]
tic()
spca <- spca(xpop, type=1, plot.nb=F, scannf=F, nfposi=20, nfnega=20)
toc()
plot(spca,1)
save(spca,file=paste0(getwd(),"/Results/Diplodus_spca_8068.RData"))

# library(interp)
# x <- other(xpop1)$xy[,1]
# y <- other(xpop1)$xy[,2]
# temp <- interp(x, y, spca$li[,1])
# image(temp, col=azur(100))

# Perform PCA
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
# load(paste0(getwd(),"/Results/Diplodus_pca.RData"))
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
    data <- data.frame(pca=spca$li[,i.axis], X=Diplodus_coord_Albers[,1], Y=Diplodus_coord_Albers[,2])
    # gstat model
    gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0))
    # Interpolation
    a <- terra::interpolate(object=pust, model=gs, xyNames=c("X", "Y"))
    pca_interp[[i.axis]] <- a[[1]]
}
# Save as raster
Diplodus_pca <- rast(pca_interp)
writeRaster(Diplodus_pca, file=paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))

# Plotting the rasters - Figure S3
g_rast <- rast(paste0(getwd(),"/Results/Diplodus_allAxes_8068.grd"))
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

# Try sPCA
load("Data_for_PCA_Mullus.RData")
tic()
spca <- spca(xpop, type=1, plot.nb=F, scannf=F, nfposi=20, nfnega=20)
toc()
# plot(spca,1)
save(spca,file=paste0(getwd(),"/Results/Mullus_spca_2753.RData"))



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
# load(paste0(getwd(),"/Results/Mullus_pca.RData"))
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
writeRaster(Mullus_pca, file=paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))

# Plotting the rasters - Figure S4
g_rast <- rast(paste0(getwd(),"/Results/Mullus_allAxes_2753.grd"))
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

