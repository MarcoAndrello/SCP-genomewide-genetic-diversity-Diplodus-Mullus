# 01b - Import conservation costs

rm(list=ls())

library(sf)
library(tidyverse)

# Load PUs
load("Planning_units.RData")

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

# Save
save(pus,file="Planning_units.RData")

