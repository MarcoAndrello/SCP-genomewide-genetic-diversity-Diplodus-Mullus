# 01 - Create planning units

rm(list=ls())

library(terra)
library(sf)
library(tidyverse)
library(units)

# Define extent of planning region (domain) from long-lat coordinates in WGS84 CRS
region <- st_sfc(
    st_polygon(
        list(rbind(c(-5.5,29), c(36.5,29), c(36.5,45), c(-5.5,45), c(-5.5,29)))
    ),
    crs=st_crs(4326))

# Transform to Albers Equal Area CRS - to get planning units of the same area (otherwise southern PUs will be larger in long-lat CRS)
region <- st_transform(region, crs=st_crs("ESRI:102013"))
# EPSG:19986 Lambers equal area recommended by the EU

# Make planning units
st_make_grid(region,cellsize=10000) %>% st_sf() -> pus
rm(region)

# Read ETOPO1 depth data
etopo1 <- rast(paste0(getwd(),"/../../Maps/ETOPO1-ice-Mediterranean.tif"))
# Extract depth from ETOPO1
pus %>% st_centroid() %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() -> pu_coord_4326
terra::extract(etopo1, pu_coord_4326) -> depth
rm(etopo1, pu_coord_4326)

# Read Fishmed data
FishMed_grid <- read.csv("FishMed_2species_1980.csv",h=T,sep=";")
FishMed_grid <- rast(FishMed_grid, type="xyz", crs="EPSG:4326")

# Extract species presence from FishMed
pus %>% st_centroid() %>% st_transform(crs=st_crs(4326)) %>% st_coordinates() %>% terra::extract(FishMed_grid, .) -> whole_Med_presence
rm(FishMed_grid)

# Add depth and species presence in grid
pus$Depth <- depth$`ETOPO1-ice-Mediterranean`
pus$Diplodus_sargus <- whole_Med_presence[,1]
pus$Mullus_surmuletus <- whole_Med_presence[,2]
rm(depth, whole_Med_presence)

# remove land cells, deep cells and cells out of domain
pus %>% filter(Depth < 0 & Depth > -200) %>% filter(!is.na(Diplodus_sargus)) -> pus
# Filtering using the column Diplodus_sargus is equivalent to filtering using the column Mullus_surmuletus:
# all(which(!is.na(pus$Diplodus_sargus)) == which(!is.na(pus$Mullus_surmuletus))) is TRUE

# Add cost, area and status
pus$cost <- 1
pus$area <- st_area(pus)
pus$status <- 0L



################################################################################
# Identify protected PUs
################################################################################

# Calculate percent protection for each PU and update status (protected or not)
mr <-
    st_read(paste0(getwd(),"/data/Abecasis_2023_MPAs/MPAFinalList.shp")) %>%
    select(NAME,Area) %>%
    st_transform(mr,crs=st_crs(pus)) 
# Total protected area
mr %>% st_area %>% sum %>% set_units("km^2")
# ## Unite overlapping MRs
mr_union <- st_union(mr)
## Add ID field
pus$ID <- 1:nrow(pus)
## Intersection: returns portions of PUs that are protected
b1 <- st_intersection(pus,st_geometry(mr_union))
## Difference: returns portions of PUs that are unprotected
b2 <- st_difference(pus,st_geometry(mr_union))
## Calculate area protected and prepare for left_join to original pus
b1$area_prot <- st_area(b1)
b1 %>% st_drop_geometry() %>% select(ID,area_prot) -> b1
## Calculate area unprotected and prepare for left_join to original pus
b2$area_unprot <- st_area(b2)
b2 %>% st_drop_geometry() %>% select(ID,area_unprot) -> b2
## Left join b1.1 and b2.1 to original pus
pus %>% left_join(b1, by="ID") %>% left_join(b2,by="ID") -> pus ###
## Replace NA areas with 0 m^2
pus$area_prot[which(is.na(pus$area_prot))] <- 0
pus$area_unprot[which(is.na(pus$area_unprot))] <- 0
## Calculate percent protection for each PU
pus$perc_prot <- pus$area_prot / pus$area
pus$perc_prot <- units::drop_units(pus$perc_prot)

## Set "status" to 1 for perc_prot > 0 or 0.01, or 0.1
pus %>% mutate(status_0 = 0, status_0.01 = 0, status_0.1 = 0) -> pus
pus$status_0[which(pus$perc_prot>0)] <- 1
pus$status_0.01[which(pus$perc_prot>0.01)] <- 1
pus$status_0.1[which(pus$perc_prot>0.1)] <- 1
rm(b1,b2,mr, mr_union)
pus %>% mutate(status_0.5 = 0) -> pus
pus$status_0.5[which(pus$perc_prot>0.5)] <- 1

# Number of protected PUs and area protected by those, according to different thresholds
pus %>% filter(status_0==1) %>% pull(area_prot) %>% length
pus %>% filter(status_0.01==1) %>% pull(area_prot) %>% length
pus %>% filter(status_0.1==1) %>% pull(area_prot) %>% length
pus %>% filter(status_0.5==1) %>% pull(area_prot) %>% length
pus %>% filter(status_0==1) %>% pull(area_prot) %>% sum %>% set_units("km^2")
pus %>% filter(status_0.01==1) %>% pull(area_prot) %>% sum %>% set_units("km^2")
pus %>% filter(status_0.1==1) %>% pull(area_prot) %>% sum %>% set_units("km^2")
pus %>% filter(status_0.5==1) %>% pull(area_prot) %>% sum %>% set_units("km^2")



################################################################################
# Import conservation costs
################################################################################

# Read shapefile and values 
grid <- st_read(paste0(getwd(),"/data/Mazor et al 2014 Opportunity costs/GridMEDSEA_Clip_etrs89.shp"))
values <- read.csv(paste0(getwd(),"/data/Mazor et al 2014 Opportunity costs/CombinedCostSenario7.csv"),sep=";")

# Join values to the sf
grid %>% 
    left_join(values,by="id") %>% 
    select(id,cost) %>% 
    st_transform(st_crs(pus)) -> 
    grid

# Intersect the two datasets and calculate the area of each intersection
st_intersection(pus, grid) %>% 
    mutate(AREA = st_area(.)) %>% 
    select(ID, AREA, cost.1) -> 
    intersection

# Calculate the mean cost as weighted means, the weights are the areas
intersection %>% 
    mutate(AREA_NUM = as.numeric(AREA)) %>% 
    group_by(ID) %>% 
    st_drop_geometry() %>% 
    summarise(mean_cost = weighted.mean(cost.1,AREA_NUM)) ->
    mean_cost

# Join the mean costs to the PUs
left_join(pus, mean_cost, by = "ID") -> pus

# Calculate coordinates of the centroid of each PU
pus %>%
    st_geometry() %>%
    st_centroid() %>%
    st_coordinates() -> pus_centroid

# Save
save(pus,pus_centroid,file="Planning_units.RData")

