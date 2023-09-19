rm(list=ls())
library(sf)
library(tidyverse)
set.seed(123)

pus <- st_read("reticolo_esempio_raptr.gpkg")
pus %>% select(-c(id, left, top, right, bottom)) -> pus
pus$Axis1 <- 1.5 + rnorm(nrow(pus), sd=2)
pus$Axis2 <- 1 + rnorm(nrow(pus), sd=2)
plot(pus)
hist(pus$Axis1)
x1 <- hist(pus$Axis1)
x1$mids
x1$density
cut(pus$Axis1,breaks=x1$breaks)
x2 <- hist(pus$Axis2)
x2$mids

# Demand points (defined by hand)
dp <- data.frame(x=c(0,2,5),y=c(1,2,3))
library(RANN)
nn2(st_drop_geometry(pus),
    dp)

sel1 <- c(9,17,6)
seg_sel1 <- data.frame(x0=dp$x, y0=dp$y,
                       x1=pus$Axis1[sel1], y1=pus$Axis2[sel1])
pus %>%
    st_drop_geometry() %>%
    plot(pch=15, col="gray")
points(st_drop_geometry(pus)[sel1,],col="green", pch=15)
points(dp,pch=16)
segments(seg_sel1$x0, seg_sel1$y0, seg_sel1$x1, seg_sel1$y1,col="blue", lwd=2)

sel2 <- c(25,2,36)
seg_sel2 <- data.frame(x0=dp$x, y0=dp$y,
                       x1=pus$Axis1[sel2], y1=pus$Axis2[sel2])
pus %>%
    st_drop_geometry() %>%
    plot(pch=15, col="gray")
points(st_drop_geometry(pus)[sel2,],col="green", pch=15)
points(dp,pch=16)
segments(seg_sel2$x0, seg_sel2$y0, seg_sel2$x1, seg_sel2$y1,col="blue", lwd=2)

plot(pus["Axis1"])
plot(pus["Axis2"])

plot(st_geometry(pus))
plot(st_geometry(pus)[sel1],add=T,border="green", lwd=3)
plot(st_geometry(pus))
plot(st_geometry(pus)[sel2],add=T,border="green", lwd=3)


