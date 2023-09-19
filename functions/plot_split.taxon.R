### PLOT THE EXAMPLE FIGURE
plot_split.taxon <- function(pus_split.taxon, Diplodus_ng_values) {
  # Read land sf
  st_read(paste0(getwd(),"/../../Maps/ne_10m_land.shp")) %>%
    st_transform(st_crs(pus_split.taxon)) %>%
    st_geometry() %>% 
    st_crop(st_bbox(pus_split.taxon)) ->
    land
  # Attach continuous sPCA values in the sf
  pus_split.taxon$Diplodus_neutr_Axis1 <- Diplodus_ng_values[,1]
  pus_split.taxon$Diplodus_neutr_Axis1[which(pus_split.taxon$Diplodus_sargus==0)] <- NA
  pus_split.taxon$which.class <- "0"
  pus_split.taxon$which.class[which(pus_split.taxon$Diplodus_neutr_Axis1_1==1)] <- "1"
  pus_split.taxon$which.class[which(pus_split.taxon$Diplodus_neutr_Axis1_2==1)] <- "2"
  pus_split.taxon$which.class[which(pus_split.taxon$Diplodus_neutr_Axis1_3==1)] <- "3"
  # Figure: continuous values
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_cont.png"),width = 10, height=5, units="cm",res=600)
  par(mar=c(1,1,1,1))
  plot(land,col="darkgray", lwd=0.5)
  plot(st_geometry(pus_split.taxon),border=NA,col="lightgray",add=T)
  plot(pus_split.taxon["Diplodus_neutr_Axis1"],border=NA,
       pal=heat.colors(100,rev=T), breaks=seq(-3,5,0.08),add=T)
  dev.off()
  # Figure: legend continuous values
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_cont_legend.png"),width = 10, height=5, units="cm",res=600)
  plot(pus_split.taxon["Diplodus_neutr_Axis1"],border=NA,
       pal=heat.colors(100,rev=T), breaks=seq(-3,5,0.08), key.pos=1)
  dev.off()
  # Figure: split taxon into three classes
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon.png"),width = 10, height=5, units="cm",res=600)
  par(mar=c(1,1,1,1))
  plot(land,col="darkgray", lwd=0.5)
  plot(pus_split.taxon["which.class"],border=NA,pal=c("lightgray","yellow","orange","red"),
       add=T)
  dev.off()
  # Figure: legend three classes
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_legend.png"),width = 10, height=5, units="cm",res=600)
  plot(pus_split.taxon["which.class"],border=NA,pal=c("lightgray","yellow","orange","red"), key.pos=1)
  dev.off()
  # Figure: split taxon: class 1
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_1.png"),width = 10, height=5, units="cm",res=600)
  par(mar=c(1,1,1,1))
  plot(land,col="darkgray", lwd=0.5)
  plot(pus_split.taxon["Diplodus_neutr_Axis1_1"],border=NA,pal=c("lightgray","yellow"),add=T)
  dev.off()
  # Figure: split taxon: class 2
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_2.png"),width = 10, height=5, units="cm",res=600)
  par(mar=c(1,1,1,1))
  plot(land,col="darkgray", lwd=0.5)
  plot(pus_split.taxon["Diplodus_neutr_Axis1_2"],border=NA,pal=c("lightgray","orange"),add=T)
  dev.off()
  # Figure: split taxon: class 3
  png(paste0(getwd(),"/Figures_split_taxon/Figure_split.taxon_3.png"),width = 10, height=5, units="cm",res=600)
  par(mar=c(1,1,1,1))
  plot(land,col="darkgray", lwd=0.5)
  plot(pus_split.taxon["Diplodus_neutr_Axis1_3"],border=NA,pal=c("lightgray","red"),add=T)
  dev.off()
  ### END PLOT
}
