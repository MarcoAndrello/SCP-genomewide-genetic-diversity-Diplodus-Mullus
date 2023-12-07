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
