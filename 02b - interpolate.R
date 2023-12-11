# PCA interpolation

rm(list=ls())

library(tidyverse)
library(terra)
library(sf)
library(gstat)
library(adegenet)

load("Planning_units.RData")
# Create raster covering the extent of pus (for interpolation)
pust <- rast(pus, resolution = c(10000, 10000), crs = crs(pus))
values(pust) <- 1
# Use only the occupied PUs
mask(pust, (pus %>% filter(Mullus_surmuletus == 1))) -> pusr


load("Data_for_PCA_Mullus.RData")
load("Results/Mullus_spca_2753.RData")
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
Diplodus_coord_Albers <- st_coordinates(obs)
rm(obs)

# Function for RMSE
RMSE <- function(observed, predicted) {
    sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

null <- perf_nn <- vector()
rmse_nn <- array(NA,dim=c(40,5))
nnmsk <- list()
for (i.axis in 1: 40) {
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=spca$li[,i.axis], X=Diplodus_coord_Albers[,1], Y=Diplodus_coord_Albers[,2])
    # 0 - Null model
    null[i.axis] <- RMSE(data$pca, mean(data$pca))
    null
    # 1 - Nearest neighbor
    gs <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp = 0))
    nn <- interpolate(pust, gs, xyNames=c("X", "Y"), debug.level=0)
    nnmsk[[i.axis]] <- mask(nn,pus)
    # k-fold cross-validation
    set.seed(20231211)
    kf <- sample(1:5, nrow(data), replace=TRUE)
    for (k in 1:5) {
        test <- data[kf == k, ]
        train <- data[kf != k, ]
        gscv <- gstat(formula=pca~1, locations=~X+Y, data=train, nmax=5, set=list(idp = 0))
        p <- predict(gscv, test, debug.level=0)$var1.pred
        rmse_nn[i.axis, k] <- RMSE(test$pca, p)
    }
    # relative model performance
    perf_nn[i.axis] <- 1 - (mean(rmse_nn[i.axis,]) / null[i.axis])
}

rmse_nn
barplot(perf_nn)
# we retain only perf > 0 (negative means the null model is better)



# 
# # 2 - Inverse distance weighting
# gs <- gstat(formula=pca~1, locations=~X+Y, data=data)
# idw <- interpolate(pust, gs, xyNames=c("X", "Y"), debug.level=0)
# idwmsk <- mask(idw, pus)
# plot(idwmsk, 1)
# points(datav,col=terrain.colors(10))
# rmse_idw <- rep(NA, 5)
# for (k in 1:5) {
#     test <- data[kf == k, ]
#     train <- data[kf != k, ]
#     gscv <- gstat(formula=pca~1, locations=~X+Y, data=train)
#     p <- predict(gscv, test, debug.level=0)$var1.pred
#     rmse_idw[k] <- RMSE(test$pca, p)
# }
# mean(rmse_idw)
# # relative model performance
#  1 - (mean(rmse_idw) / null)
# round(perf_idw, 3)
# 
# # Fit a variogram
# v <- variogram(gs,width=50000)
# v
# plot(v)
# 
# # Optimization
# f1 <- function(x, test, train) {
#     nmx <- x[1]
#     idp <- x[2]
#     if (nmx < 1) return(Inf)
#     if (idp < .001) return(Inf)
#     m <- gstat(formula=pca~1, locations=~X+Y, data=train, nmax=nmx, set=list(idp=idp))
#     p <- predict(m, newdata=test, debug.level=0)$var1.pred
#     RMSE(test$pca, p)
# }
# set.seed(20150518)
# i <- sample(nrow(data), 0.2 * nrow(data))
# tst <- data[i,]
# trn <- data[-i,]
# opt <- optim(c(8, .5), f1, test=tst, train=trn)
# str(opt)
# 
# # Best model
# m <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=opt$par[1], set=list(idp=opt$par[2]))
# idw <- interpolate(pust, m, xyNames=c("X", "Y"), debug.level=0)
# idwmsk <- mask(idw, pus)
# plot(idw,1)
# plot(idwmsk, 1)
# points(datav,col=terrain.colors(10))
# values(idwmsk) %>% as.vector %>% unique %>% length
# 
# 
# # 3 - Inverse path distance weighting
# library(ipdw)
# SeaRaster <- rast("C:/Users/marco/Desktop/seamap05.tif")
# plot(SeaRaster)
# plot(project(datav, crs(SeaRaster)),add=T)
# # Generate a stack of path accumulated distance raster objects
# datav %>% st_as_sf %>% 
#     # .[1:2,] %>%
#     st_transform(crs(SeaRaster))  -> data2
# # SeaRaster %>% project(data2) %>% as("Raster") -> sea_rast
# rstack <- pathdistGen(data2, as(SeaRaster,"Raster"), 100, progressbar = T)
# plot(rstack,20)
# plot(data2[20,],add=T)
# 
# final.raster <- ipdwInterp(data2, rstack, paramlist = c("pca"), 
#                            dist_power = 0.5, overlapped = TRUE)
# plot(final.raster)
# plot(data2, add = TRUE)
# plot(data2)
# 
# 
