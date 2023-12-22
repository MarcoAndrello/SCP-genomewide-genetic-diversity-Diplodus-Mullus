# Interpolation study

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


####################
# Diplodus sargus
####################
# Load xpop and results of the PCA
load("Data_for_PCA_Diplodus.RData")
load("Results/Diplodus_pca_pop.RData")

# Interpolation
# Define obs sf object so we can convert the lat-lon coordinates to Albers equal area
obs <- st_as_sf(data.frame(pca1=0,
                           lon=xpop@other$xy$x,
                           lat=xpop@other$xy$y),
                coords=c("lon","lat"),
                crs=4326)
obs <- st_transform(obs,crs=st_crs("ESRI:102013"))
Diplodus_coord_Albers <- st_coordinates(obs)
rm(obs)

# Object to store the Mean Absolute Prediction Error for all axes, the 5 folds and the two interpolation methods
mape_nn <- mape_idw <- array(NA,dim=c(17,5))
set.seed(20231214)
kf <- sample(1:5, nrow(Diplodus_coord_Albers), replace=TRUE)
# Loop on PCA axes
for (i.axis in 1 : 17){ #ncol(pca$li)){
    if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=pca$li[,i.axis], X=Diplodus_coord_Albers[,1], Y=Diplodus_coord_Albers[,2])
    # k-fold for nn
    for (k in 1:5) {
        test <- data[kf == k, ]
        train <- data[kf != k, ]
        gscv <- gstat(formula=pca~1, locations=~X+Y, data=train, nmax=5, set=list(idp=0))
        p <- predict(gscv, test)$var1.pred
        mape_nn[i.axis,k] <- mean(abs((test$pca-p)/test$pca)) * 100
    }
    # k-fold for idw
    for (k in 1:5) {
        test <- data[kf == k, ]
        train <- data[kf != k, ]
        gscv <- gstat(formula=pca~1, locations=~X+Y, data=train)
        p <- predict(gscv, test)$var1.pred
        mape_idw[i.axis,k] <- mean(abs((test$pca-p)/test$pca)) * 100
    }
}

# Plot boxplots showing the distribution of MAPE over the 5 folds, for each Axis and for the two methods (for Diplodus)
colnames(mape_nn) <- colnames(mape_idw) <- c(1:5) 
mape_nn %>% as_tibble %>% mutate(method="nn") %>%               # MAPE of nn
    rbind(mape_idw %>% as_tibble %>% mutate(method="idw")) %>%  # MAPE of idw
    mutate(Axis=rep(c(1:17),2)) %>%
    pivot_longer(cols=c(1:5),names_to="kfold") %>%
    rename(MAPE=value) %>%
    ggplot(aes(x=factor(Axis),y=MAPE)) +
    geom_boxplot(aes(fill=method))

# Plot map of nn interpolation for Axis 1
data <- data.frame(pca=pca$li[,1], X=Diplodus_coord_Albers[,1], Y=Diplodus_coord_Albers[,2])
gsnn <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0))
ann <- terra::interpolate(object=pust, model=gsnn, xyNames=c("X", "Y"))
annmsk <- mask(ann,filter(pus, Diplodus_sargus==1))
plot(annmsk,1,main="Nearest neighbor interpolation")
values(annmsk) %>% as.vector %>% unique %>% length

# Plot map of idw interpolation for Axis 1
gsidw <- gstat(formula=pca~1, locations=~X+Y, data=data)
aidw <- terra::interpolate(object=pust, model=gsidw, xyNames=c("X", "Y"))
aidwmsk <- mask(aidw,filter(pus, Diplodus_sargus==1))
plot(aidwmsk,1, main="Inverse distance weighting interpolation")
values(aidwmsk) %>% as.vector %>% unique %>% length

## NOTA: RIFARE  IL GRAFICO E LE DUE MAPPE PRECEDENTI SALVANDOLE, E METTENDOLE IN MATERIALI SUPPLEMENTARI
## COMMENTARE IL CODICE SOTTO PER MULLUS


####################
# Mullus surmuletus
####################

load("Data_for_PCA_Mullus.RData")
load("Results/Mullus_pca_pop.RData")

# Interpolation
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
i.axis <- 1
mape_nn <- mape_idw <- array(NA,dim=c(26,5))
set.seed(20231214)
kf <- sample(1:5, nrow(coord_Albers), replace=TRUE)
for (i.axis in 1 : 26){
    if (i.axis %% 10 == 0) cat(i.axis,"\n"); flush.console()
    # Dataframe with observation (pca scores) and X and Y coord
    data <- data.frame(pca=pca$li[,i.axis], X=coord_Albers[,1], Y=coord_Albers[,2])
    # k-fold nn
    for (k in 1:5) {
        test <- data[kf == k, ]
        train <- data[kf != k, ]
        gscv <- gstat(formula=pca~1, locations=~X+Y, data=train, nmax=5, set=list(idp=0))
        p <- predict(gscv, test)$var1.pred
        mape_nn[i.axis,k] <- mean(abs((test$pca-p)/test$pca)) * 100
    }
    # k-fold idw
    for (k in 1:5) {
        test <- data[kf == k, ]
        train <- data[kf != k, ]
        gscv <- gstat(formula=pca~1, locations=~X+Y, data=train)
        p <- predict(gscv, test)$var1.pred
        mape_idw[i.axis,k] <- mean(abs((test$pca-p)/test$pca)) * 100
    }
}

mape_nn %>% rowMeans %>% barplot(xlab="PCA axis",names.arg=c(1:26), main="MAPE Nearest neighbor",ylim=c(0,15000))
mape_idw %>% rowMeans %>% barplot(xlab="PCA axis",names.arg=c(1:26), main="MAPE inverse distance weighing",ylim=c(0,15000))

colnames(mape_nn) <- colnames(mape_idw) <- c(1:5) 
mape_nn %>% as_tibble %>% mutate(method="nn") %>%
    rbind(mape_idw %>% as_tibble %>% mutate(method="idw")) %>%
    mutate(Axis=rep(c(1:26),2)) %>%
    pivot_longer(cols=c(1:5),names_to="kfold") %>%
    rename(MAPE=value) %>%
    ggplot(aes(x=factor(Axis),y=MAPE)) +
    geom_boxplot(aes(fill=method)) +
    ylim(0,2500)

#
data <- data.frame(pca=pca$li[,1], X=coord_Albers[,1], Y=coord_Albers[,2])
gsnn <- gstat(formula=pca~1, locations=~X+Y, data=data, nmax=5, set=list(idp=0))
ann <- terra::interpolate(object=pust, model=gsnn, xyNames=c("X", "Y"))
annmsk <- mask(ann,filter(pus, Mullus_surmuletus==1))
plot(annmsk,1,main="Nearest neighbor interpolation")
values(annmsk) %>% as.vector %>% unique %>% length

gsidw <- gstat(formula=pca~1, locations=~X+Y, data=data)
aidw <- terra::interpolate(object=pust, model=gsidw, xyNames=c("X", "Y"))
aidwmsk <- mask(aidw,filter(pus, Mullus_surmuletus==1))
plot(aidwmsk,1, main="Inverse distance weighting interpolation")
values(aidwmsk) %>% as.vector %>% unique %>% length

