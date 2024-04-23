### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, terra, and sf packages for R.
###
### source('C:/Ecology/Drive/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/01_preparing_data.r')
### source('E:/Adam/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/01_preparing_data.r')
###
### CONTENTS ###
### setup ###
### obtain shapes for countries with focal basins ###
### obtain elevation data from AWS Terrain Tiles ###

#############
### setup ###
#############

	rm(list=ls())
	
	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	setwd(paste0(drive, '/Research/fasterRaster - Streamlined GIS in R through GRASS'))
	
	# library(fasterRaster) # GIS
	# devtools::load_all(paste0(drive, '/R/fasterRaster')) # GIS
	library(elevatr) # elevation
	library(enmSdmX) # GIS/species distribution modeling
	library(geodata) # GIS data
	library(omnibus) # utilities
	library(terra) # rasters and vectors

	source('./fasterRaster_benchmarking/00_constants.r')

	faster(grassDir = 'C:/Program Files/GRASS GIS 8.3', memory = 1024 * 8)

# say('#####################################################')
# say('### obtain shapes for countries with focal basins ###')
# say('#####################################################')

# 	gadm <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/gadm_410.gpkg'))
# 	seAsia <- gadm[gadm$NAME_0 %in% focal_country_names]

# 	for (country in focal_country_names) {

# 		say(country)

# 		x <- seAsia[seAsia$NAME_0 == country]
# 		x <- aggregate(x)

# 		path <- paste0('C:/Ecology/Research Data/GADM/Version 4.1/', country, '.gpkg')
# 		writeVector(x, path, overwrite = TRUE)

# 	}

# say('####################################################')
# say('### obtain elevation data from AWS Terrain Tiles ###')
# say('####################################################')

# 	library(elevatr)

# 	se_asian_basins <- vect('./data/Hydrological Basins in Southeast Asia (FAO)/hydrobasins_asia.shp')
# 	focal_basins <- se_asian_basins[se_asian_basins$MAJ_NAME %in% priary_basins_large]
# 	focal_extent <- ext(focal_basins)
# 	focal_extent <- as.polygons(focal_extent, crs = crs(focal_basins))
# 	focal_extent <- buffer(focal_extent, study_region_buffer_size_m)
# 	focal_extent <- st_as_sf(focal_extent)

# 	# elevation from AWS Terrain Tiles, zoom 11 (~32 m)
# 	elev_aws_z11_m <- get_elev_raster(focal_extent, z = 11, override_size_check = TRUE)
# 	elev_aws_z11_m <- round(elev_aws_z11_m, 0L)
# 	elev_aws_z11_m <- setMinMax(elev_aws_z11_m)
# 	names(elev_aws_z11_m) <- 'elev_aws_z11_m'
# 	writeRaster(elev_aws_z11_m, './data/AWS Terrain Tiles/elev_aws_z11_m.tif', overwrite = TRUE, datatype = 'INT2S')

# 	# elevation from AWS Terrain Tiles, zoom 7 (~512 m)
# 	elev_aws_z7_m <- get_elev_raster(focal_extent, z = 7, override_size_check = TRUE)
# 	elev_aws_z7_m <- round(elev_aws_z7_m, 0L)
# 	elev_aws_z7_m <- setMinMax(elev_aws_z7_m)
# 	names(elev_aws_z7_m) <- 'elev_aws_z7_m'
# 	writeRaster(elev_aws_z7_m, './data/AWS Terrain Tiles/elev_aws_z7_m.tif', overwrite = TRUE, datatype = 'INT2S')
