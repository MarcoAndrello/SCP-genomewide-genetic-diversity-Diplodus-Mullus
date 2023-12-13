# 02 Perform spatial PCA and interpolates the results

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(adegenet)
library(gstat)
library(tictoc)
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

# Function for RMSE
RMSE <- function(observed, predicted) {
    sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

for (i.species in 1 : 2) {
    
    name_species = c("Diplodus", "Mullus")[i.species]
    
    if (name_species == "Diplodus") {
        gendata_filename = "data/Diplodus_8068.RData"
        coord_filename = "Diplodus_coord.RData"
        res_readdata_filename = "Data_for_sPCA_Diplodus.RData"
        res_spca_filename = "Results/Diplodus_spca_8068.RData"
        colname = "Diplodus_sargus"
        res_interp_filename = "Results/Diplodus_allAxes_8068.grd"
        png_filename_1 = "Genetic_rasters_Diplodus.png"
    }
    if (name_species == "Mullus") {
        gendata_filename = "data/Mullus_2753.RData"
        coord_filename = "Mullus_coord.RData"
        res_readdata_filename = "Data_for_sPCA_Mullus.RData"
        res_spca_filename = "Results/Mullus_spca_2753.RData"
        colname = "Mullus_surmuletus"
        res_interp_filename = "Results/Mullus_allAxes_2753.grd"
        png_filename_1 = "Genetic_rasters_Mullus_1.png"
        png_filename_2 = "Genetic_rasters_Mullus_2.png"
    }
    
    # Load merged (neutral + adaptive) dataset (but this excludes the non-outlier that are not in HWE)
    load(gendata_filename)
    if (name_species == "Diplodus") {
        x <- Diplodus_8068; rm(Diplodus_8068)
        x$pop <- Diplodus_sampling$SamplingCell 
    } else {
        x <- Mullus_2753; rm(Mullus_2753)
        # Add population information
        x$pop <- Mullus_sampling$SamplingCell 
    }
    # Convert to genpop
    xpop <- genind2genpop(x)
    # Match rownames of xpop to sampling cells in cell_sampling to retrieve geographic coordinates
    match_xpop_to_cellsampling <- match(row.names(xpop$tab),cell_sampling$SamplingCell)
    # Retrieve geographic coordinates
    coord <- cell_sampling[match_xpop_to_cellsampling,c("Longitude","Latitude")]
    rm(match_xpop_to_cellsampling)
    save(coord,file=coord_filename)
    # Add spatial coordinates
    other(xpop)$xy <- coord 
    names(other(xpop)$xy) <- c("x","y")
    save(xpop,file=res_readdata_filename)
    
    # Perform sPCA: 1.5 hrs for Diplodus. 2 min for Mullus
    load(res_readdata_filename)
    tic()
    spca <- spca(xpop, type=1, plot.nb=F, scannf=F, nfposi=20, nfnega=20)
    toc()
    save(spca,file=res_spca_filename)
    screeplot(spca)# (takes several minutes to plot)
    
    load(res_spca_filename)
    ## Transform coordinates from 4326 to Albers projection
    obs <- st_as_sf(data.frame(pca1=0,
                               lon=xpop@other$xy$x,
                               lat=xpop@other$xy$y),
                    coords=c("lon","lat"),
                    crs=4326)
    obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
    coord_Albers <- st_coordinates(obs)
    rm(obs)
    
    # Spatial interpolation: loop on axes
    null <- perf_nn <- vector()
    rmse_nn <- array(NA,dim=c(40,5))
    nnmsk <- list()
    for (i.axis in 1 : 40) {
        # Dataframe with observation (spca scores) and X and Y coord
        data <- data.frame(spca=spca$li[,i.axis], X=coord_Albers[,1], Y=coord_Albers[,2])
        # 0 - Null model
        null[i.axis] <- RMSE(data$spca, mean(data$spca))
        null
        # 1 - Nearest neighbor
        # k-fold cross-validation
        set.seed(20231211)
        kf <- sample(1:5, nrow(data), replace=TRUE)
        for (k in 1:5) {
            test <- data[kf == k, ]
            train <- data[kf != k, ]
            gscv <- gstat(formula=spca~1, locations=~X+Y, data=train, nmax=5, set=list(idp = 0))
            p <- predict(gscv, test, debug.level=0)$var1.pred
            rmse_nn[i.axis, k] <- RMSE(test$spca, p)
        }
        # relative model performance
        perf_nn[i.axis] <- 1 - (mean(rmse_nn[i.axis,]) / null[i.axis])
        # Interpolate on the full domain
        if (perf_nn[i.axis] < 0) next
        gs <- gstat(formula=spca~1, locations=~X+Y, data=data, nmax=5, set=list(idp = 0))
        nn <- interpolate(pust, gs, xyNames=c("X", "Y"), debug.level=0)
        # Mask the domain to areas of species presence
        nnmsk[[i.axis]] <- mask(nn,filter(pus, (!!sym(colname))==1))$var1.pred
    }
    summary(perf_nn[which(perf_nn>0)])
    spca_rast <- rast(nnmsk[which(perf_nn>0)])
    writeRaster(spca_rast, file=res_interp_filename)
    
    # Plotting the rasters - Figure S4
    # spca_rast <- rast(res_interp_filename)
    names(spca_rast) <- paste("Axis",c(1:nlyr(spca_rast))) # not really the "original" axes, but those filtered out by RMSE criterion
    # Read countries for plotting maps
    ne_countries(scale = 50, returnclass = "sf") %>% st_transform(st_crs(spca_rast)) -> countries
    png(png_filename_1,width=20,height=21,units="cm",res=300)
    par(mfrow=c(5,3))
    for (i in 1 : min(nlyr(spca_rast),15)) {
        plot(spca_rast,i,mar=c(0.5,0.1,1,0.1),axes=F,legend=F,main=paste("Axis",i))
        polys(countries,col="gray",lwd=0.01)
    }
    dev.off()
    if (i.species == 1) break
    png(png_filename_2,width=20,height=21,units="cm",res=300)
    par(mfrow=c(5,3))
    for (i in 16 : nlyr(spca_rast)) {
        plot(spca_rast,i,mar=c(0.5,0.1,1,0.1),axes=F,legend=F,main=paste("Axis",i))
        polys(countries,col="gray",lwd=0.01)
    }
    dev.off()
    
}

