rm(list=ls())

library(tidyverse)
library(sf)
library(wdpar)

# Planning units
load("Planning_units.RData")

v.p <- c("Pelagos","Santuario per i Mammiferi Marini","Tutela del Tursiops truncatus")

iso_codes <- c("ESP",
               "FRA",
               "MCO",
               "ITA",
               "SVN",
               "HRV",
               "BIH",
               "MNE",
               "ALB",
               "GRC",
               "TUR",
               "SYR",
               "LBN",
               "ISR",
               "PSE",
               "EGY",
               "LBY",
               "TUN",
               "DZA",
               "MAR",
               "MLT",
               "CYP")

mpas <- list()

# pas <- wdpa_fetch("ITA", download_dir=paste0(getwd(),"/ita-pa"),check_version = F)

for (i.country in 15 : length(iso_codes)){
    cat(iso_codes[i.country],"\n")
    pas <- wdpa_fetch(iso_codes[i.country], download_dir=paste0(getwd(),"/existing-pa"),wait=T)
    
    # Keep only coastal and marine
    pas %>% filter(MARINE > 0) -> mpas_raw
    
    # for ITA, FRA and Monaco: remove Pelagos sanctuary
    for (i.p in 1 : length(v.p)){
        if (length(grep(v.p[i.p],mpas_raw$NAME))>0) mpas_raw <-mpas_raw[-grep(v.p[i.p],mpas_raw$NAME),]
        # mpas_raw <-mpas_raw[-grep("Pelagos",mpas_raw$NAME),]
        # mpas_raw <-mpas_raw[-grep("Santuario per i Mammiferi Marini",mpas_raw$NAME),]
        # mpas_raw <-mpas_raw[-grep("Tutela del Tursiops truncatus",mpas_raw$NAME),]
    }
    
    if (nrow(mpas_raw) == 0) next
    ## 
    # Cleaning
    mpas[[i.country]] <- wdpa_clean(mpas_raw, crs = st_crs(pus)$input, erase_overlaps = F)
    # rm(mpas_raw)
    # st_write(mpas,dsn="mpas.gpkg")
    
    # # Remove MPAs less than 1km2 surface area
    # mpas[-which(mpas$AREA_KM2<1),] %>% nrow()
}

# rbind mpas...
mpas_binder <- bind_rows(mpas)

save(mpas_binder,file="mpas_binder.RData")

load("mpas_binder.RData")
mpas <- mpas_binder
rm(mpas_binder)

st_write(mpas,dsn="MPAs_2099.gpkg")

# Further filtering
mpas %>% 
    filter(DESIG_ENG != "Special Protection Area (Birds Directive)") %>%
    filter(DESIG_ENG != "Special Protection Areas") %>%
    filter(DESIG_ENG != "Specialy Protected Area") %>%
    filter(NAME != "Corredor de migración de cetáceos del Mediterráneo") %>%
    filter(NAME != "Grands dauphins de l'Agriate") -> mpas

st_write(mpas,dsn="MPAs_1702.gpkg")

#  STOP HERE
# THE FOLLOWING IS TO CREATE A DIFFERENT PUS FILE WITH SEPARATE POLYGONS FOR PROTECTED AND UNPROTECTED PORTIONS OF THE SAME PU

# Dissolve to unite overlapping MPAs
mpasd <- wdpa_dissolve(mpas)

# Intersection with PUs:
# results in a sf with PU within MPAs and PU outside MPAs. no more square PUs.
# plot(st_geometry(pus))
# plot(st_geometry(mpasd),add=T)
# Intersection: returns portions of PUs that are protected
b1 <- st_intersection(pus,st_geometry(mpasd))
# Difference: returns portions of PUs that are unprotected
b2 <- st_difference(pus,st_geometry(mpasd))
# is.pa argument say which geometries are protected
b1$is.pa <- T
b2$is.pa <- F
# plot(st_geometry(mpasd))
# plot(st_geometry(b1),add=T,border="blue")
# plot(st_geometry(b2),add=T,border="red")
pus_MPA <- rbind(b1,b2)
pus_MPA$area <- st_area(pus_MPA)
nrow(pus_MPA)
nrow(pus_MPA[pus_MPA$is.pa==T,])
plot(pus_MPA["is.pa"],border=NA)
sum(pus_MPA$area[pus_MPA$is.pa==T]) *1e-6
sum(pus_MPA$area[pus_MPA$is.pa==T]) / sum(pus_MPA$area)

save(pus_MPA,file="Planning_units_MPA.RData")
st_write(pus_MPA, dsn="Planning_units_MPA.gpkg",append=F)

# Remove small holes: try debuffer-buffer 100m
pus_MPA$AREA_KM2 <- st_area(pus_MPA) *1e-6
summary(pus_MPA$AREA_KM2)
pus_MPA_debuf <- st_buffer(pus_MPA,dist= -100)
pus_MPA_debuf$AREA_KM2 <- st_area(pus_MPA_debuf) * 1e-6
summary(pus_MPA_debuf$AREA_KM2)
pus_MPA_buf <- st_buffer(pus_MPA_debuf,dist= 100)
pus_MPA_buf$AREA_KM2 <- st_area(pus_MPA_buf) * 1e-6
summary(pus_MPA_buf$AREA_KM2)
which(drop_units(pus_MPA_buf$AREA_KM2)==0)
