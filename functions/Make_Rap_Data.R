# Function to make rap_data for my scenario analysis
# Now works only for Diplodus sargus

Make_Rap_Data <- function(pus,
                          species = "Diplodus_sargus",
                          make.geographic.space = F,
                          make.neutral.space = F,
                          make.outlier.space = F,
                          n_demand.point_geog = 50,
                          # n_demand.point_neutral = 50L,   # Not used bcs we place demand points in each cell of the regular grid covering the space
                          # n_demand.point_outlier = 50L,
                          # quantile_demand.point_neutral = 0.75, # Not used bcs we place demand points in each cell of the regular grid covering the space
                          # quantile_demand.point_outlier = 0.75,
                          ngridcells_demand.point_neutral = 10,
                          ngridcells_demand.point_outlier = 10,
                          ng = NULL,
                          og = NULL,
                          target_amount = 0.3,
                          target_geog = 0.3,
                          target_neutral = 0.99,
                          target_outlier = 0.99,
                          plot=F
) {
    
    # Create species probabilities. For 2 species: rbind(...,...)
    pu.species.probabilities <- data.frame(species=1L,
                                           pu=c(1:nrow(pus)),
                                           value=pus[[species]])
    pu.species.probabilities <- filter(pu.species.probabilities, value==1)
    
    # Coonvert to PolySet
    pus %>%  sf::as_Spatial() %>% raptr::SpatialPolygons2PolySet() -> polygons
    # pus %>%  convert2PolySet() -> polygons
    
    # Calculate boundary
    boundary <- raptr::calcBoundaryData(polygons)
    
    
    # Attribute Spaces for the geographic space
    if (make.geographic.space == T) {
        cat("Geographic space not implemented\n")
        make.geographic.space <- F
        # ## keep only polygon of presence, unite them (dissolve)
        # pus %>% filter(Diplodus_sargus == 1) %>% st_geometry() %>% st_union() -> FishMed_pol_Diplodus
        # ## Sample n_demand.point_geog points within area of presence (thanks to Jeff Hanson for this idea)
        # dp_geographic_sf <- st_sample(FishMed_pol_Diplodus,n_demand.point_geog)
        # ## Plot to see if it's OK
        # # plot(FishMed_pol_Diplodus)
        # # plot(dp_geographic_sf,add=T,col="red",pch=16)
        # ## Define geographic attribute space (it will be a list with 2 elements when having 2 species)
        # geographic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=st_coordinates(st_centroid(st_geometry(pus))), ids=1:nrow(pus)),
        #                                          demand.points = DemandPoints(coords=st_coordinates(dp_geographic_sf), weights=rep(1,length(dp_geographic_sf))),
        #                                          species=1L))
        # rm(FishMed_pol_Diplodus, dp_geographic_sf)
    }
    
    # Attribute Spaces for the neutral genetic space
    if (make.neutral.space == T) {
        
        # Keeping only the PUs occupied by the species
        pus %>%
            filter(.data[[species]] == 1) %>%
            ## Extracting coordinates of PUs in neutral genetic space
            st_geometry() %>%
            st_centroid() %>%
            st_coordinates() %>%
            raster::extract(x=ng,y=.) ->
            pu_coord_neutr.gen.space 
        ## Convert to sf object
        pu_coord_neutr.gen.space %>%
            as.data.frame() %>%
            st_as_sf(coords=c(1,2)) -> pu_coord_neutr.gen.space_points
        ## Make a grid covering the extent of the Coordinates of PUs in this space
        pu_coord_neutr.gen.space_points %>%
            st_make_grid(n=c(ngridcells_demand.point_neutral, ngridcells_demand.point_neutral)) -> grid
        # Calculate how many PUs fall into each grid cell
        grid %>%
            st_intersects(pu_coord_neutr.gen.space_points) %>%
            lengths() -> numpoints
        ## Convert the grid to sf
        grid_sf <- st_sf(NUMPOINTS=numpoints, geometry=grid)
        ## Keep only non-empty grid cells
        grid_sf %>%
            filter(NUMPOINTS>0) ->
            grid_sf_temp
        cat("Number of demand points for neutral genetic space:",nrow(grid_sf_temp),"\n")
        ## Transform grid to points
        grid_sf_temp %>%
            st_geometry() %>%
            st_centroid() ->
            st_geometry(grid_sf_temp)
        ## Set the coordinates of points as demand points, with weight = number of planning units
        dp_neutral <- DemandPoints(coords = st_coordinates(grid_sf_temp),
                                   weights = grid_sf_temp$NUMPOINTS/max(grid_sf_temp$NUMPOINTS))
        
        ## Plot
        if (plot==T) plot_demand_points_neutral(grid_sf, pu_coord_neutr.gen.space_points, dp_neutral) 
        
        ## Define neutral genetic space (it will be a list with 2 elements when having 2 species)
        neutral_genetic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=pu_coord_neutr.gen.space, ids=pu.species.probabilities$pu),
                                                      demand.points = dp_neutral,
                                                      species=1L))
    }
    
    
    # Attribute Spaces for the outlier genetic space
    if (make.outlier.space == T) {
        
        # Keeping only the PUs occupied by the species
        pus %>%
            filter(.data[[species]] == 1) %>%
            ## Extracting coordinates of PUs in outlier genetic space
            st_geometry() %>%
            st_centroid() %>%
            st_coordinates() %>%
            raster::extract(x=og,y=.) ->
            pu_coord_outlier.gen.space 
        ## Convert to sf object
        pu_coord_outlier.gen.space %>%
            as.data.frame() %>%
            st_as_sf(coords=c(1,2)) -> pu_coord_outlier.gen.space_points
        ## Make a grid covering the extent of the Coordinates of PUs in this space
        pu_coord_outlier.gen.space_points %>%
            st_make_grid(n=c(ngridcells_demand.point_outlier, ngridcells_demand.point_outlier)) ->
            grid
        # Calculate how many PUs fall into each grid cell
        grid %>%
            st_intersects(pu_coord_outlier.gen.space_points) %>%
            lengths() -> numpoints
        ## Convert the grid to sf
        grid_sf <- st_sf(NUMPOINTS=numpoints, geometry=grid)
        ## Keep only non-empty grid cells
        grid_sf %>%
            filter(NUMPOINTS>0) ->
            grid_sf_temp
        cat("Number of demand points for outlier genetic space:",nrow(grid_sf_temp),"\n")
        ## Transform grid to points
        grid_sf_temp %>%
            st_geometry() %>%
            st_centroid() ->
            st_geometry(grid_sf_temp)
        ## Set the coordinates of points as demand points, with weight = number of planning units
        dp_outlier <- DemandPoints(coords = st_coordinates(grid_sf_temp),
                                   weights = grid_sf_temp$NUMPOINTS/max(grid_sf_temp$NUMPOINTS))
        ## Plot
        if (plot==T) plot_demand_points_outlier(grid_sf, pu_coord_outlier.gen.space_points, dp_outlier) 
        
        ## Define outlier genetic space (it will be a list with 2 elements when having 2 species)
        outlier_genetic_spaces <- list(AttributeSpace(planning.unit.points = PlanningUnitPoints(coords=pu_coord_outlier.gen.space, ids=pu.species.probabilities$pu),
                                                      demand.points = dp_outlier,
                                                      species=1L))
        
    }
    
    # Define attribute.spaces
    attribute.spaces <- list()
    if (make.geographic.space == T) {
        attribute.spaces[[(length(attribute.spaces)+1)]] <-
            AttributeSpaces(
                spaces = geographic_spaces,
                name = "geographic")
    }
    if (make.neutral.space == T) {
        attribute.spaces[[(length(attribute.spaces)+1)]] <-
            AttributeSpaces(
                spaces = neutral_genetic_spaces,
                name = "neutral_genetic")
    }
    if (make.outlier.space == T) {
        attribute.spaces[[(length(attribute.spaces)+1)]] <-
            AttributeSpaces(
                spaces = outlier_genetic_spaces,
                name = "outlier_genetic")
    }
    
    
    # Create targets
    ## Geographic targets to 0
    targets <- data.frame(species = 1L,
                          target = 0:3L,
                          proportion = c(target_amount,
                                         target_geog,
                                         target_neutral,
                                         target_outlier),
                          included = c(TRUE,
                                       make.geographic.space,
                                       make.neutral.space,
                                       make.outlier.space))
    
    targets %>% dplyr::filter(included) %>% dplyr::select(-included) -> targets
    targets$target <- seq(0,(nrow(targets)-1))
    
    rap_data <- RapData(pu = st_drop_geometry(pus),
                        species = data.frame(name=species),
                        targets = targets,
                        pu.species.probabilities = pu.species.probabilities,
                        attribute.spaces = attribute.spaces,
                        boundary = boundary,
                        polygons = polygons)
    return(rap_data)
}