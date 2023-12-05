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

