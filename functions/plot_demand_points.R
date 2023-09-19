
plot_demand_points_neutral <- function(grid_sf, pu_coord_neutr.gen.space_points, dp_neutral) {
    
    # Figure a
    png("Figure Xa.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="white",reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2")
    plot(pu_coord_neutr.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure b
    png("Figure Xb.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2",
         col="white")
    plot(pu_coord_neutr.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure c
    fillcolors <- heat.colors(10,rev=T)[cut(grid_sf$NUMPOINTS,breaks=seq(0,500,50),include.lowest=F)]
    png("Figure Xc.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2",
         col=fillcolors)
    plot(pu_coord_neutr.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure c for legend
    png("Figure Xc legend.png", width=15, height=15, units="cm", res=600)
    plot(grid_sf, border="gray",reset=F, axes=T,
         pal=heat.colors(10,rev=T),breaks=seq(0,500,50),key.pos=1)
    dev.off()
    
    # Figure d
    dp_neutral_sf <- st_as_sf(data.frame(dp_neutral@coords, weights=dp_neutral@weights),
                              coords = c("X","Y"))
    png("Figure Xd.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2",
         col=fillcolors)
    plot(pu_coord_neutr.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    plot(st_geometry(dp_neutral_sf), reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2",
         pch=16, col="red", cex=dp_neutral_sf$weights*2.5,
         add=T)
    dev.off()
    
    # Figure d for legend
    png("Figure Xd legend.png", width=15, height=15, units="cm", res=600)
    plot(rep(1,5), pch=16, col="red", cex=seq(0.2,1,0.2)*2.5,xlim=c(-1,+7))
    dev.off()
    
}


plot_demand_points_outlier <- function(grid_sf, pu_coord_outlier.gen.space_points, dp_outlier) {
    
    # Figure a
    png("Figure Ya.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="white",reset=F, axes=T,
         xlab = "Outlier sPCA Axis 1", ylab="Outlier sPCA Axis 2")
    plot(pu_coord_outlier.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure b
    png("Figure Yb.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Outlier sPCA Axis 1", ylab="Outlier sPCA Axis 2",
         col="white")
    plot(pu_coord_outlier.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure c
    fillcolors <- heat.colors(10,rev=T)[cut(grid_sf$NUMPOINTS,breaks=seq(0,500,50),include.lowest=F)]
    png("Figure Yc.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Outlier sPCA Axis 1", ylab="Outlier sPCA Axis 2",
         col=fillcolors)
    plot(pu_coord_outlier.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    dev.off()
    
    # Figure c for legend
    png("Figure Yc legend.png", width=15, height=15, units="cm", res=600)
    plot(grid_sf, border="gray",reset=F, axes=T,
         pal=heat.colors(10,rev=T),breaks=seq(0,500,50),key.pos=1)
    dev.off()
    
    # Figure d
    dp_outlier_sf <- st_as_sf(data.frame(dp_outlier@coords, weights=dp_outlier@weights),
                              coords = c("X","Y"))
    png("Figure Yd.png", width=15, height=15, units="cm", res=600)
    plot(st_geometry(grid_sf), border="gray",reset=F, axes=T,
         xlab = "Outlier sPCA Axis 1", ylab="Outlier sPCA Axis 2",
         col=fillcolors)
    plot(pu_coord_outlier.gen.space_points,
         pch=21, bg=scales::alpha("black",0.1), col="black",
         add=T)
    plot(st_geometry(dp_outlier_sf), reset=F, axes=T,
         xlab = "Neutral sPCA Axis 1", ylab="Neutral sPCA Axis 2",
         pch=16, col="red", cex=dp_outlier_sf$weights*2.5,
         add=T)
    dev.off()
    
    # Figure d for legend
    png("Figure Yd legend.png", width=15, height=15, units="cm", res=600)
    plot(rep(1,5), pch=16, col="red", cex=seq(0.2,1,0.2)*2.5,xlim=c(-1,+7))
    dev.off()
    
}

