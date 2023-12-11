# VARIO CODICE PER INTERPOLAZIONI SPAZIALI

#########################################################################
# RSPAT E GSPATIAL
#########################################################################


# https://rspatial.org/analysis/4-interpolation.html

remotes::install_github('rspatial/rspat')

library(rspat)
d <- spat_data('precipitation')
head(d)

mnts <- toupper(month.abb)
d$prec <- rowSums(d[, mnts])
plot(sort(d$prec), ylab="Annual precipitation (mm)", las=1, xlab="Stations")

dsp <- vect(d, c("LONG", "LAT"), crs="+proj=longlat +datum=NAD83")
CA <- spat_data("counties")
# define groups for mapping
cuts <- c(0,200,300,500,1000,3000)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'darkblue'))
plot(CA, col="lightgray", lwd=4, border="darkgray")
plot(dsp, "prec", type="interval", col=blues(10), legend=TRUE, cex=2,
     breaks=cuts, add=TRUE, plg=list(x=-117.27, y=41.54))
lines(CA)

TA <- "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=WGS84 +units=m"
dta <- project(dsp, TA)
cata <- project(CA, TA)

# Null model
RMSE <- function(observed, predicted) {
    sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
null <- RMSE(dsp$prec, mean(dsp$prec))
null

# proximity polygons
v <- voronoi(dta)
plot(v)
points(dta)

vca <- crop(v, cata)
plot(vca, "prec")
points(dta, col=blues(10), breaks=cuts)

r <- rast(vca, res=10000)
vr <- rasterize(vca, r, "prec")
plot(vr)

set.seed(5132015)
kf <- sample(1:5, nrow(dta), replace=TRUE)
rmse <- rep(NA, 5)
for (k in 1:5) {
    test <- dta[kf == k, ]
    train <- dta[kf != k, ]
    v <- voronoi(train)
    p <- extract(v, test)
    rmse[k] <- RMSE(test$prec, p$prec)
}
rmse
## [1] 192.0568 203.1304 183.5556 177.5523 205.6921
mean(rmse)
## [1] 192.3974
# relative model performance
perf <- 1 - (mean(rmse) / null)
round(perf, 3)
## [1] 0.558


library(gstat)
d <- data.frame(geom(dta)[,c("x", "y")], as.data.frame(dta))
head(d)
gs <- gstat(formula=prec~1, locations=~x+y, data=d, nmax=5, set=list(idp = 0))
nn <- interpolate(r, gs, debug.level=0)
nnmsk <- mask(nn, vr)
plot(nnmsk, 1)

rmsenn <- rep(NA, 5)
for (k in 1:5) {
    test <- d[kf == k, ]
    train <- d[kf != k, ]
    gscv <- gstat(formula=prec~1, locations=~x+y, data=train, nmax=5, set=list(idp = 0))
    p <- predict(gscv, test, debug.level=0)$var1.pred
    rmsenn[k] <- RMSE(test$prec, p)
}
rmsenn
mean(rmsenn)
1 - (mean(rmsenn) / null)

#########################################################################
# spatialRF Moran Eigenvector Maps
#########################################################################

# Try MEMs
library(spatialRF)
library(units)
datav %>% st_as_sf %>% st_geometry %>% st_distance %>% as.matrix %>% drop_units -> dist_mat
mem(dist_mat)
data$dummy <- 1
m1 <- rf_spatial(data=data, dependent.variable.name = "pca", predictor.variable.name = "dummy",
                distance.matrix = dist_mat, method="mem.moran.sequential")
## Now with in-sea distance. Calculate the distance
library(gdistance)
# Create transition matrices based on conductance rasters
# Load the raster of marine areas
# This should be already a conductance layer: Coastal Sea = 1, others (land and abyss) = 0
SeaRaster <- raster::raster("C:/Users/marco/Desktop/seamap.tif")
# Create transition matrix and correct it for lat-long
t <- transition(SeaRaster, function(x) mean(x), directions=8)
t <- geoCorrection(t, type="c")
    
# Calculate least coast distances
datav %>% st_as_sf %>% st_geometry %>% st_transform(st_crs(4326)) %>% st_coordinates -> fromCoords
sea_dist <-  costDistance(t, fromCoords) 
### guardare anche terra::costDist
sea_dist %>% as.matrix -> sea_dist
m2 <- rf_spatial(data=data, dependent.variable.name = "pca", predictor.variable.name = "dummy",
                distance.matrix = sea_dist, method="mem.moran.sequential")



#########################################################################
# inverse path distance weighting
#########################################################################


rm(list=ls())

library(tidyverse)
library(sf)
library(ipdw)

pols <- st_read(system.file("extdata/kattegat_coast.gpkg", package = "ipdw"))
pnts <- st_read(system.file("extdata/kattegat_pnts.gpkg", package = "ipdw"))

costras <- costrasterGen(pnts, pols, extent = "pnts",
                         projstr = projection(pols))
# insert contiguous barrier
costras[160:170, 1:80] <- 10000

library(spatstat)

W              <- owin(range(c(st_bbox(pnts)["xmin"], st_bbox(pnts)["xmax"])),
                       range(c(st_bbox(pnts)["ymin"], st_bbox(pnts)["ymax"])))
kat.pp         <- ppp(st_coordinates(pnts)[,1], st_coordinates(pnts)[,2], window = W)
mean.neighdist <- mean(nndist(kat.pp))

# grid building
gridsize       <- mean.neighdist * 2
grainscale.fac <- gridsize / res(costras)[1]
gridras        <- aggregate(costras, fact = grainscale.fac)
gridpol        <- rasterToPolygons(gridras)
gridpol$value  <- row.names(gridpol)

# spatial join
fulldataset.over <- sf::st_join(pnts, st_as_sf(gridpol))

# grid selection
set.seed(2)
gridlev <- unique(fulldataset.over$value)
for (i in seq_along(gridlev)) {
    activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
    selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
    if (i == 1) {
        training <- activesub[selectnum, ]
    } else {
        training <- rbind(training, activesub[selectnum, ])
    }
}


validate <- fulldataset.over[!(row.names(fulldataset.over) %in%
                                   row.names(training)), ]
plot(costras)
plot(st_geometry(training), add = TRUE)
plot(st_geometry(validate), col = "red", add = TRUE)


paramlist <- c("salinity")
final.ipdw <- ipdw(training, costras, range = mean.neighdist * 10, paramlist,
                   overlapped = TRUE)
final.ipdw2 <- ipdw(training, costras, range = mean.neighdist*20, paramlist, dist_power = 0.5,
                    overlapped = TRUE)

plot(final.ipdw, main = "Kattegat salinity (ppt)")
plot(final.ipdw2, main = "Kattegat salinity (ppt) 2")
plot(pnts,add=T)

##### ottenere il raster delle distances
library(sf)
sf_ob <- data.frame(rnorm(2))
xy    <- data.frame(x = c(4, 2), y = c(8, 4))
sf_ob <- st_as_sf(cbind(sf_ob, xy), coords = c("x", "y"))

m <- matrix(NA, 10, 10)
costras <- raster(m, xmn = 0, xmx = ncol(m), ymn = 0, ymx = nrow(m))

# introduce spatial gradient
costras[] <- runif(ncell(costras), min = 1, max = 10)
for (i in 1:nrow(costras)) {
    costras[i, ] <- costras[i, ] + i
    costras[, i] <- costras[, i] + i
}

rstack <- pathdistGen(sf_ob, costras, 100, progressbar = FALSE)



#########################################################################
# phylin
#########################################################################
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