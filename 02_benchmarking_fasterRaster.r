### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, fasterRaster, and sf packages for R.
###
### source('C:/Ecology/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/02_benchmarking_fasterRaster.r')
### source('E:/Adam/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/02_fasterRaster_benchmarking_fasterRaster.r')
###
### CONTENTS ###
### settings ###
### benchmark fasterRaster ###

################
### settings ###
################

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	setwd(paste0(drive, '/Research/fasterRaster - Streamlined GIS in R through GRASS'))

	library(enmSdmX) # SDMing
	library(omnibus) # utilities
	library(ROCR) # AUC
	library(speedglm) # GLM
	library(terra) # GLM

	fff <- function() devtools::load_all(paste0(drive, '/R/fasterRaster'))
	fff()
	source(paste0(drive, '/R/fasterRaster/R/rbind.r'))
	# library(fasterRaster) # rasters and vectors

	library(usdm) # predictor selection

	# memory <- 16 * 1024
	memory <- 56 * 1024
	# cores <- 2
	cores <- 4

	work_dir <- paste0('C:/!scratch/grass')
	grass_dir <- 'C:/Program Files/GRASS GIS 8.3/'
	faster(grassDir = grass_dir, memory = memory, cores = cores, workDir = work_dir)

say('##############################')
say('### benchmark fasterRaster ###')
say('##############################')

	### demesne
	###########

	# demesne <- 'Large'
	# basins_primary <- primary_basins_large
	# basins_secondary <- NA
	# country_names <- focal_country_names

	# demesne <- 'Medium'
	# basins_primary <- 'Salween'
	# basins_secondary <- NA
	# country_names <- c('China', 'Myanmar', 'Thailand')

	demesne <- 'Small'
	basins_primary <- 'Mekong'
	basins_secondary <- 'Nam Loi'
	country_names <- c('China', 'Myanmar')

	### start
	#########

	sink(paste0('./outputs/', tolower(demesne), '_fasterRaster_benchmark.txt'), split = TRUE)
	say('BENCHMARK FASTERRASTER')
	say(date(), post = 1)

	source('./fasterRaster_benchmarking/00_constants.r')

	say('SESSION INFO', level = 2)
	print(sessionInfo())

	say('FASTERRASTER OPTIONS', level = 2)
	print(faster())

	say('DEMESNE', level = 2)
	say('Demesne: ....................................................... ', demesne)
	say('Primary Basins: ................................................ ', paste(basins_primary, collapse = ', '))
	say('Secondary Basins: .............................................. ', paste(basins_secondary, collapse = ', '))
	say('Country Names: ................................................. ', paste(country_names, collapse = ', '), post = 1)

	say('SETTINGS', level = 2)
	say('Drive with datasets and libraries: ............................. ', drive)
	say('Number of cores used for multi-core functions: ................. ', cores)
	say('Maximum memory allowed (GB): ................................... ', memory / 1024)
	say('Study region buffer size (m): .................................. ', study_region_buffer_size_m)
	say('Number of folds: ............................................... ', n_folds)
	say('Proportion of forest loss cells to sample per fold: ............ ', cross_valid_prop)
	say('Quantiles used to divide final prediction raster into zones: ... ', paste(threshold_quantiles, collapse = ', '))
	say('Number of iterations in permutation importance test: ........... ', n_permute)
	say('')

	all_ok <- function() {
		x <- readline('Are the settings above OK? (y/n) ')
		if (x != 'y') stop('All your base are belong to us.')
	}
	
	all_ok()

	# stores runtime information
	timings <- data.frame()

	### LOAD RAW RASTERS
	####################
	say('LOAD RAW RASTERS', level = 2)
	step <- 'Load raw rasters'

	# forest cover for 2000 and 2020
	tile_names <- if (demesne == 'Large') {
		c('30N_090E', '20N_100E', '20N_090E', '10N_100E', '40N_090E', '30N_100E')
	} else if (demesne == 'Medium') {
		c('30N_090E', '40N_090E', '20N_090E', '30N_100E')
	} else if (demesne == 'Small') {
		c('30N_090E', '30N_100E')
	}
	forest_files_2000 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2020 Forest Height & Extent/2000/', tile_names, '.tif')
	forest_files_2020 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2020 Forest Height & Extent/2020/', tile_names, '.tif')
	
	forest_2000 <- list()
	start <- Sys.time()
	for (i in seq_along(forest_files_2000)) {
		forest_2000[[i]] <- fast(forest_files_2000[i])
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = '2000 forest cover', dtype = 'raster', start = start, stop = stop, n = length(forest_2000))
	
	forest_2020 <- list()
	start <- Sys.time()
	for (i in seq_along(forest_files_2020)) {
		forest_2020[[i]] <- fast(forest_files_2020[i])
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = '2020 forest cover', dtype = 'raster', start = start, stop = stop, n = length(forest_2020))
	
	# elevation rasters
	start <- Sys.time()
	elev_scale_11_m <- fast('./data/AWS Terrain Tiles/elev_aws_z11_m.tif')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	elev_scale_7_m <- fast('./data/AWS Terrain Tiles/elev_aws_z7_m.tif')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Elevation at scale 7', dtype = 'raster', start = start, stop = stop)
	
	# LULC
	lulc_files_2000 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2022 LULC/2000/', tile_names, '.tif')
	lulc_levels <- read.csv(paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2022 LULC/legend.csv'))

	lulc <- list()
	start <- Sys.time()
	for (i in seq_along(lulc_files_2000)) {
		lulc[[i]] <- fast(lulc_files_2000[i])
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = length(lulc))
	
	start <- Sys.time()
	for (i in seq_along(lulc_files_2000)) {
		levels(lulc[[i]]) <- lulc_levels
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'levels()', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = length(lulc))

	# human population
	start <- Sys.time()
	population <- fast(paste0(drive, '/Research Data/Population - Gridded Population of the World (GPW) v4/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_30_sec.tif'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	# soil nitrogen
	start <- Sys.time()
	nitrogen <- fast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean.tif'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Soil nitrogen', dtype = 'raster', start = start, stop = stop)

	# soil organic carbon
	start <- Sys.time()
	soc <- fast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean.tif'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	### LOAD RAW VECTORS
	####################
	say('LOAD RAW VECTORS', level = 2)
	step <- 'Load raw vectors'

	# drainage basins
	start <- Sys.time()
	se_asia_basins <- fast(paste0(drive, '/Research Data/Watershed Basins in Southeast Asia (FAO)/hydrobasins_asia.gpkg'), verbose = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Watershed basins', dtype = 'vector', start = start, stop = stop)

	# countries of SE Asia
	countries <- list()
	start <- Sys.time()
	for (i in seq_along(country_names)) {

		country <- country_names[i]
		countries[[i]] <- fast(paste0(drive, '/Research Data/GADM/Version 4.1/', country, '.gpkg'), verbose = TRUE)

	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'SE Asian countries', dtype = 'vector', start = start, stop = stop, n = length(countries))

	# rivers
	start <- Sys.time()
	rivers <- fast(paste0(drive, '/Research Data/Rivers - FAO/rivers_asia_37331.shp'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Asian rivers', dtype = 'vector', start = start, stop = stop)
	
	# roads vectors
	roads <- list()
	start <- Sys.time()
	for (i in seq_along(country_names)) {
	
		country <- country_names[i]
		say('roads: ', country)
		if (country == 'India') country <- paste0(country, '-north-eastern-zone')
		country <- tolower(country)
		fn <- paste0(drive, '/Research Data/Open Street Maps (Geofabrik)/', country, '-latest-free.shp/gis_osm_roads_free_1.shp')
		roads[[i]] <- fast(fn, dropTable = TRUE, snap = 0, area = 0, verbose = TRUE)
	
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Roads', dtype = 'vector', start = start, stop = stop, n = length(roads))
	
	### DEFINE STUDY REGION
	#######################
	say('DEFINE STUDY REGION', level = 2)
	step <- 'Define study region'
	
	### select study region
	start <- Sys.time()
	if (demesne == "Small") {
		focal_basins <- se_asia_basins[se_asia_basins$MAJ_NAME %in% basins_primary & se_asia_basins$SUB_NAME %in% basins_secondary]
	} else {
		focal_basins <- se_asia_basins[se_asia_basins$MAJ_NAME %in% basins_primary]
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Basins vector', dtype = 'vector', start = start, stop = stop)
	
	### aggregate sub-basins (need to add small buffer first to ensure overlap)
	start <- Sys.time()
	focal_basins_buff <- buffer(focal_basins, width = 0.01) # buffer in degrees
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'buffer()', target = 'Mekong basin vector', dtype = 'vector', start = start, stop = stop)
	
	start <- Sys.time()
	focal_basins_agg <- aggregate(focal_basins_buff, dissolve = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'aggregate()', target = 'Basins vector', dtype = 'vector', start = start, stop = stop)

	# unbuffer basins to create study region
	start <- Sys.time()
	study_region <- buffer(focal_basins_agg, width = -0.01)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'buffer()', target = 'Aggregated basins vector', dtype = 'vector', start = start, stop = stop)

	### buffer study region
	# need to project to equal-area first to get buffer in meters, unlike in terra::buffer()
	start <- Sys.time()
	study_region_equal_area <- project(study_region, getCRS('Asia Lambert'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	start <- Sys.time()
	study_region_equal_area_buffered <- buffer(study_region_equal_area, study_region_buffer_size_m)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'buffer()', target = 'Basins vector', dtype = 'vector', start = start, stop = stop)

	# project back to WGS84
	start <- Sys.time()
	study_region_buffered <- project(study_region_equal_area_buffered, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	### FOREST COVER
	################
	say('FOREST COVER', level = 2)
	step <- 'Prepare forest cover rasters'

	# merge
	start <- Sys.time()
	forest_2000 <- do.call(merge, forest_2000)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'merge()', target = '2000 forest cover', dtype = 'raster', start = start, stop = stop, n = length(forest_2000))

	start <- Sys.time()
	forest_2020 <- do.call(merge, forest_2020)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'merge()', target = '2020 forest cover', dtype = 'raster', start = start, stop = stop, n = length(forest_2020))
	
	# rasterize study region
	# NB We are doing this to use as a mask. We can also use mask(), but multiplying by this raster is faster than repeated calls to mask().
	start <- Sys.time()
	study_region_rast <- rasterize(study_region_buffered, forest_2000)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rasterize()', target = 'Study region and 2000 forest cover', dtype = 'raster/vector', start = start, stop = stop)

	# trim study region
	start <- Sys.time()
	study_region_rast <- trim(study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'trim()', target = 'Study region raster', dtype = 'raster', start = start, stop = stop)

	# aggregated study region
	# aggregate raster to make processing faster
	start <- Sys.time()
	study_region_rast_agg <- aggregate(study_region_rast, fact = 64)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'aggregate()', target = 'Study region raster', dtype = 'raster', start = start, stop = stop)
	
	# crop
	start <- Sys.time()
	forest_2000 <- crop(forest_2000, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = '2000 forest cover & study region', dtype = 'raster/vector', start = start, stop = stop)
	
	start <- Sys.time()
	forest_2020 <- crop(forest_2020, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = '2020 forest cover & study region', dtype = 'raster/vector', start = start, stop = stop)
	
	# mask
	start <- Sys.time()
	forest_2000 <- forest_2000 * study_region_rast
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = '2000 forest cover & study region raster', dtype = 'raster', start = start, stop = stop)
	
	start <- Sys.time()
	forest_2020 <- forest_2020 * study_region_rast
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = '2020 forest cover & study region raster', dtype = 'raster', start = start, stop = stop)
	
	# name (do not need a name for 2020 forest)
	start <- Sys.time()
	names(forest_2000) <- 'forest_2000'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = '2000 forest cover', dtype = 'raster', start = start, stop = stop)

	# forest change
	start <- Sys.time()
	forest_delta <- forest_2000 - forest_2020
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = '2000 and 2020 forest cover', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	forest_loss <- forest_delta == 1
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'comparison', target = 'Forest change raster', dtype = 'raster', start = start, stop = stop)
	
	# mask forest loss raster so it contains only cells that had forest in 2000
	start <- Sys.time()
	forest_2000_no_forest_masked <- app(forest_2000, fun = '= if (forest_2000 == 0, null(), forest_2000)')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'app()', target = 'Forest 2000 raster', dtype = 'raster', start = start, stop = stop)
	
	start <- Sys.time()
	forest_loss <- forest_loss * forest_2000_no_forest_masked
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Masked forest 2000 & loss raster', dtype = 'raster', start = start, stop = stop)
	
	# name
	start <- Sys.time()
	names(forest_loss) <- 'forest_loss'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)

	### FOREST DENSITY
	##################
	say('FOREST DENSITY', level = 2)
	step <- 'Prepare forest density predictors'

	### calculate forest density @ scale 1
	start <- Sys.time()
	forest_density_3cells <- focal(forest_2000, w = 3, fun = 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'focal()', target = '2000 forest cover at 3-cell scale', dtype = 'raster', start = start, stop = stop)

	# calculate forest density @ scale 2
	start <- Sys.time()
	forest_density_33cells <- focal(forest_2000, w = 33, fun = 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'focal()', target = '2000 forest cover at 33-cell scale', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(forest_density_3cells) <- 'forest_density_3cells'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Forest density at 3-cell scale', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	names(forest_density_33cells) <- 'forest_density_33cells'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Forest density at 33-cell scale', dtype = 'raster', start = start, stop = stop)

	### FOREST FRAGMENTATION CLASS
	##############################
	say('FOREST FRAGMENTATION CLASS', level = 2)
	step <- 'Prepare forest fragmentation class predictor'

	start <- Sys.time()
	forest_frag_class <- fragmentation(forest_2000)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fragmentation()', target = '2000 forest cover', dtype = 'raster', start = start, stop = stop)

	# names
	start <- Sys.time()
	names(forest_frag_class) <- 'forest_frag_class'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Fragmentation raster', dtype = 'raster', start = start, stop = stop)
	
	### ELEVATION
	#############
	say('ELEVATION', level = 2)
	step <- 'Prepare elevation predictor'

	# resample
	start <- Sys.time()
	elev_scale_7_m <- resample(elev_scale_7_m, forest_loss, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Elevation at scale 7 & forest loss', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	elev_scale_11_m <- resample(elev_scale_11_m, forest_loss, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Elevation at scale 11 & forest loss', dtype = 'raster', start = start, stop = stop)

	# mask with study region
	start <- Sys.time()
	elev_scale_7_m <- mask(elev_scale_7_m, study_region_buffered)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Elevation at scale 7 & study region', dtype = 'raster/vector', start = start, stop = stop)

	start <- Sys.time()
	elev_scale_11_m <- mask(elev_scale_11_m, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Elevation at scale 11 & study region raster', dtype = 'raster', start = start, stop = stop)

	# names
	start <- Sys.time()
	names(elev_scale_11_m) <- 'elev_scale_11_m'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	### SLOPE
	#########
	say('SLOPE', level = 2)
	step <- 'Prepare slope predictor'

	# slope
	start <- Sys.time()
	slope_scale_7 <- terrain(elev_scale_7_m, 'slope')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'terrain()', target = 'Elevation at scale 7', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	slope_scale_11 <- terrain(elev_scale_11_m, 'slope')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'terrain()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(slope_scale_7) <- 'slope_scale_7'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Slope at scale 7', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	names(slope_scale_11) <- 'slope_scale_11'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Slope at scale 11', dtype = 'raster', start = start, stop = stop)

	### DISTANCE TO RIVERS
	######################
	say('DISTANCE TO RIVERS', level = 2)
	step <- 'Prepare distance to rivers predictor'

	# crop
	start <- Sys.time()
	rivers <- crop(rivers, study_region_buffered)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Rivers', dtype = 'vector', start = start, stop = stop)
	
	# distance
	start <- Sys.time()
	dist_to_rivers_agg_km <- distance(study_region_rast_agg, rivers, unit = 'km')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'distance()', target = '2000 forest cover (coarse scale) & rivers', dtype = 'raster/vector', start = start, stop = stop)
	
	# resample
	start <- Sys.time()
	dist_to_rivers_km <- resample(dist_to_rivers_agg_km, forest_2000, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = '2000 forest cover (coarse and fine scale)', dtype = 'raster', start = start, stop = stop)

	# mask
	start <- Sys.time()
	dist_to_rivers_km <- mask(dist_to_rivers_km, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Distance to rivers', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(dist_to_rivers_km) <- 'dist_to_rivers_km'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Distance to rivers', dtype = 'raster/vector', start = start, stop = stop)
	
	### COUNTRIES
	#############
	say('COUNTRIES', level = 2)
	step <- 'Prepare country predictor'

	# crop
	start <- Sys.time()
	for (i in seq_along(countries)) {
		countries[[i]] <- crop(countries[[i]], study_region_buffered)
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Countries vectors & study region', dtype = 'vector', start = start, stop = stop, n = length(countries))

	# merge
	start <- Sys.time()
	countries <- do.call(rbind, countries)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rbind()', target = 'Countries vectors', dtype = 'vector', start = start, stop = stop, n = length(focal_country_names))

	# assign country names (carried over as factor)
	start <- Sys.time()
	countries$country <- country_names
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'dollar_assign', target = 'Countries vector', dtype = 'vector', start = start, stop = stop)

	# rasterize
	start <- Sys.time()
	countries_rast <- rasterize(countries, study_region_rast, byGeom = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rasterize()', target = 'Countries vector & study region raster', dtype = 'raster/vector', start = start, stop = stop)
	
	# name
	start <- Sys.time()
	names(countries_rast) <- 'country'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Countries raster', dtype = 'raster', start = start, stop = stop)

	### SPARSE VEGETATION
	#####################
	say('SPARSE VEGETATION', level = 2)
	step <- 'Prepare sparse vegetation predictor'

	# merge LULC rasters
	start <- Sys.time()
	lulc <- do.call(merge, lulc)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'merge()', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = length(tile_names))

	# crop
	start <- Sys.time()
	lulc <- crop(lulc, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'LULC & study region raster', dtype = 'raster', start = start, stop = stop)

	# mask
	start <- Sys.time()
	lulc <- mask(lulc, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'LULC & study region raster', dtype = 'raster', start = start, stop = stop)
	
	# mask short vegetation
	start <- Sys.time()
	short_veg <- (lulc >= 19 & lulc <= 24) | (lulc >= 102 & lulc <= 124)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'comparison', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = 4)

	# # density of short vegetation @ scale 1
	# start <- Sys.time()
	# short_veg_density_3cells <- focal(short_veg, w = 3, fun = 'sum', na.rm = TRUE)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'focal()', target = 'Short vegetation', dtype = 'raster', start = start, stop = stop)
	
	# density of short vegetation @ scale 2
	start <- Sys.time()
	short_veg_density_33cells <- focal(short_veg, w = 33, fun = 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'focal()', target = 'Short vegetation', dtype = 'raster', start = start, stop = stop)
	
	# # name
	# start <- Sys.time()
	# names(short_veg_density_3cells) <- 'short_veg_density_3cells'
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'names()', target = 'Short vegetation density at 3-cell scale', dtype = 'raster', start = start, stop = stop)
	
	start <- Sys.time()
	names(short_veg_density_33cells) <- 'short_veg_density_33cells'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Short vegetation density at 33-cell scale', dtype = 'raster', start = start, stop = stop)
	
	### ROADS
	#########
	say('ROADS', level = 2)
	step <- 'Prepare distance to roads predictor'

	# crop roads vectors
	start <- Sys.time()
	for (i in seq_along(roads)) {
		roads[[i]] <- crop(roads[[i]], study_region_buffered)
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Roads & study region', dtype = 'raster', start = start, stop = stop, n = length(roads))
	
	# combine roads vectors
	start <- Sys.time()
	roads <- do.call(rbind, roads)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rbind()', target = 'Roads', dtype = 'raster', start = start, stop = stop, n = length(focal_country_names))
	
	# distance to nearest road
	start <- Sys.time()
	dist_to_roads_agg_km <- distance(study_region_rast_agg, roads, unit = 'km')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'distance()', target = 'aggregated study region raster & roads', dtype = 'raster/vector', start = start, stop = stop)
	
	# resample
	start <- Sys.time()
	dist_to_roads_km <- resample(dist_to_roads_agg_km, forest_2000, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'distance to roads and 2000 forest raster', dtype = 'raster', start = start, stop = stop)

	# mask
	start <- Sys.time()
	dist_to_roads_km <- mask(dist_to_roads_km, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(dist_to_roads_km) <- 'dist_to_roads_km'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)
	
	### SOIL
	########
	say('SOIL', level = 2)
	step <- 'Prepare soil predictors'

	# project study region to soil CRS
	start <- Sys.time()
	study_region_igh <- project(study_region_buffered, nitrogen)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Study region & nitrogen raster', dtype = 'raster/vector', start = start, stop = stop)

	# crop soil rasters to study region
	start <- Sys.time()
	nitrogen <- crop(nitrogen, study_region_igh)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Study region & nitrogen raster', dtype = 'raster/vector', start = start, stop = stop)

	start <- Sys.time()
	soc <- crop(soc, study_region_igh)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Study region & soil organic carbon', dtype = 'raster/vector', start = start, stop = stop)

	# project soil rasters to study region CRS and forest resolution
	start <- Sys.time()
	nitrogen <- project(nitrogen, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Nitrogen & study region raster', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	soc <- project(soc, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Soil organic carbon & study region raster', dtype = 'raster', start = start, stop = stop)

	# names
	start <- Sys.time()
	names(nitrogen) <- 'nitrogen'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Nitrogen', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	names(soc) <- 'soc'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	### HUMAN POPULATION
	####################
	say('HUMAN POPULATION', level = 2)
	step <- 'Prepare human population density predictor'

	# crop
	start <- Sys.time()
	population <- crop(population, study_region_buffered)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Population density & study region', dtype = 'raster/vector', start = start, stop = stop)
	
	# resample
	start <- Sys.time()
	population <- resample(population, forest_2000, method = 'near')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Population density & forest cover in 2000', dtype = 'raster', start = start, stop = stop)

	# mask
	start <- Sys.time()
	population <- mask(population, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Population density & study region raster', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(population) <- 'population'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Population density at 3-cell scale', dtype = 'raster', start = start, stop = stop)

	### COLLATE PREDICTOR RASTERS
	#############################
	say('COLLATE PREDICTOR RASTERS', level = 2)
	step <- 'Collate predictor rasters'

	### stack continuous predictors
	start <- Sys.time()
	predictors_continuous <- c(
		# forest_density_3cells,
		forest_density_33cells,
		elev_scale_11_m,
		slope_scale_7,
		slope_scale_11,
		dist_to_rivers_km,
		dist_to_roads_km,
		# short_veg_density_3cells,
		short_veg_density_33cells,
		nitrogen,
		soc,
		population
	)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_continuous))

	# crop and mask to non-buffered study region
	start <- Sys.time()
	predictors_continuous <- crop(predictors_continuous, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_continuous))

	start <- Sys.time()
	predictors_continuous <- mask(predictors_continuous, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_continuous))

	# scale predictors
	start <- Sys.time()
	predictors_continuous <- scale(predictors_continuous)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'scale()', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_continuous))

	### non-continuous predictor rasters plus response raster
	start <- Sys.time()
	predictors_noncontinuous <- c(
		forest_loss,
		forest_frag_class,
		countries_rast
	)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_noncontinuous))

	# crop and mask to non-buffered study region
	start <- Sys.time()
	predictors_noncontinuous <- crop(predictors_noncontinuous, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictors_noncontinuous))

	### combine all continuous, non-continuous, and response rasters
	start <- Sys.time()
	env <- c(
		predictors_noncontinuous,
		predictors_continuous
	)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	# subset nitrogen and forest loss rasters to serve as mask
	# The soil rasters mask out cities, so we will use them to ensure randomly located points used to extract values for the model do not fall into NA cells.
	start <- Sys.time()
	nitrogen <- env[['nitrogen']]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_bracket', target = 'Nitrogen raster', dtype = 'raster', start = start, stop = stop)
	
	# mask to areas where there were forest cover in 2000 and soil values
	start <- Sys.time()
	soil_mask <- mask(nitrogen, focal_basins)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Nitrogen & focal basins vector', dtype = 'raster', start = start, stop = stop)
	
	# start <- Sys.time()
	# soil_mask <- soil_mask * 0 + 1
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Soil mask raster', dtype = 'raster', start = start, stop = stop, n = 2)

	start <- Sys.time()
	env <- mask(env, soil_mask)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	# retain forest loss raster for computation later
	start <- Sys.time()
	forest_loss_focal_basin <- env[['forest_loss']]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_bracket', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)
	
	# retain fragmentation class raster for computation later
	start <- Sys.time()
	forest_frag_class_focal_basin <- env[['forest_frag_class']]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_bracket', target = 'Forest fragmentation raster', dtype = 'raster', start = start, stop = stop)
	
	# rasters from which to select random points for calibration/evaluation
	# fragmentation raster has NAs around border
	start <- Sys.time()
	select_from <- c(forest_loss_focal_basin, forest_frag_class_focal_basin)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c', target = 'Forest loss & fragmentation rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(select_from))

	# save raster stack
	start <- Sys.time()
	writeRaster(env, paste0('./outputs/', tolower(demesne), '_fasterRaster_response_and_predictors.tif'), overwrite = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	### CROSS-VALIDATION SETUP
	##########################
	say('CROSS-VALIDATION SETUP', level = 2)
	step <- 'Cross-validation setup'

	# calculate number of sites for calibration and evaluation
	# we want a balanced design (equal number of sites in "no loss" and "loss" cells)
	# calculate sample size... selecting more than needed bc sites can be placed in NA cells

	start <- Sys.time()
	n_cells_loss <- global(forest_loss_focal_basin, 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'global()', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)

	n_cells_loss <- n_cells_loss$sum

	start <- Sys.time()
	n_cells_in_study_region_rast <- ncell(forest_loss_focal_basin)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'ncell()', target = 'Study region raster', dtype = 'raster', start = start, stop = stop)

	# number of random points to draw
	# inflating by 1.2 to ensure we have enough points on loss cells
	# inflating by 2 because we want a set for calibration and for evaluation
	n_points_inflated <- ceiling(2 * 1.2 * cross_valid_prop * n_cells_in_study_region_rast)

	# NB We do not select random points from the full predictor stack first. Rather, we select random points from the forest loss raster and nitrogen raster, then use it to remove points that are outside the study region (NA cells) and to select the desired number of cells of each class (forest loss/persistence). *Then* we extract from the full response/predictor stack. This is much faster than selecting points and extracting values using the full stack initially.

	# store results from modeling in this data frame
	results <- data.frame()

	for (k in 1:n_folds) {

		say('CROSS-VALIDATION FOLD ', k, level = 2)
		step <- paste0('Cross-validation: Select and extract from calibration and evaluation points')

		### randomly select points
		start <- Sys.time()
		sample_sites <- spatSample(select_from, n_points_inflated, values = TRUE, as.points = TRUE, cats = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Forest loss raster masked to focal basin', dtype = 'raster', start = start, stop = stop)

		# remove points in NA cells
		start <- Sys.time()
		values_at_points <- as.data.frame(sample_sites)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'as.data.frame()', target = 'Candidate calibration/evaluation sites', dtype = 'vector', start = start, stop = stop)
		
		completes <- which(complete.cases(values_at_points))
		start <- Sys.time()
		sample_sites <- sample_sites[completes]
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)

		# retain number of points desired
		start <- Sys.time()
		forest_loss_extract <- sample_sites$forest_loss
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_dollar', target = 'Candidate calibration/evaluation sites', dtype = 'vector', start = start, stop = stop)

		no_losses_index <- which(forest_loss_extract == 0)
		losses_index <- which(forest_loss_extract == 1)

		n_keeps <- min(round(2 * cross_valid_prop * n_cells_loss), length(no_losses_index), length(losses_index))

		if (length(no_losses_index) > n_keeps) no_losses_index <- sample(no_losses_index, n_keeps)
		if (length(losses_index) > n_keeps) losses_index <- sample(losses_index, n_keeps)

		# catch cases where number of forest loss sites != number of persistence sites
		stopifnot(length(no_losses_index) == length(no_losses_index))
		keeps <- sort(c(no_losses_index, losses_index))
		
		start <- Sys.time()
		sample_sites <- sample_sites[keeps]
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)

		### extract from full response/predictor stack
		start <- Sys.time()
		response_predictors <- extract(env, sample_sites, cats = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'extract()', target = 'Predictor and response rasters & calibration/evaluation sites', dtype = 'raster/vector', start = start, stop = stop, n = length(env))

		# catch cases with NAs
		if (any(!complete.cases(response_predictors))) stop('At least one cell has NAs.')

		# only conduct variable selection on first fold... use same variables for remainder of folds
		if (k == 1) {

			### select variables with low collinearity
			# predictors_continuous_df <- response_predictors[ , c('forest_density_3cells', 'forest_density_33cells', 'elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'dist_to_rivers_km', 'dist_to_roads_km', 'short_veg_density_3cells', 'short_veg_density_33cells', 'nitrogen', 'soc', 'population')]
			predictors_continuous_df <- response_predictors[ , c('forest_density_33cells', 'elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'dist_to_rivers_km', 'dist_to_roads_km', 'short_veg_density_33cells', 'nitrogen', 'soc', 'population')]

			n <- nrow(predictors_continuous_df)
			cors <- vifcor(predictors_continuous_df, method = 'spearman', size = n, th = 0.7)

			say('Collinearity analysis and variable selection:', pre = 1)
			print(cors)
			say('')
			
			continuous_predictors_selected <- cors@results$Variables

			vars_selected <- c('forest_loss', continuous_predictors_selected, 'forest_frag_class', 'country')
			predictors_selected <- c(continuous_predictors_selected, 'forest_frag_class', 'country')

		} # if first fold

		# add factor predictors
		response_predictors <- response_predictors[ , ..vars_selected]

		# divide into calibration/evaluation
		# Split loss/no loss into their own data frames, then assign random folds, then equilibrate so we have an equal number of cases in each fold in each loss/no loss group.
		response_predictors$fold <- NA_integer_		
		response_predictors_loss <- response_predictors[response_predictors$forest_loss == 1, ]
		response_predictors_no_loss <- response_predictors[response_predictors$forest_loss == 0, ]
		
		# apportion folds equally across LOSS sites
		n <- nrow(response_predictors_loss)
		fold <- sample(0L:1L, n, replace = TRUE)
		ones <- which(fold == 1)
		zeroes <- which(fold == 0)
		
		n_ones <- length(ones)
		n_zeroes <- length(zeroes)

		if (n_ones > n_zeroes) {
		
			diff <- n_ones - n_zeroes
			half_diff <- (diff / 2)
			swaps <- sample(ones, half_diff)
			ones <- ones[-which(ones %in% swaps)]
			zeroes <- c(zeroes, swaps)
		
		} else if (n_ones < n_zeroes) {
		
			diff <- n_zeroes - n_ones
			half_diff <- (diff / 2)
			swaps <- sample(zeroes, half_diff)
			zeroes <- zeroes[-which(zeroes %in% swaps)]
			ones <- c(ones, swaps)
		
		}
		response_predictors_loss$fold[ones] <- 1
		response_predictors_loss$fold[zeroes] <- 0

		# apportion folds equally across NO LOSS sites
		n <- nrow(response_predictors_no_loss)
		fold <- sample(0L:1L, n, replace = TRUE)
		ones <- which(fold == 1)
		zeroes <- which(fold == 0)
		
		n_ones <- length(ones)
		n_zeroes <- length(zeroes)

		if (n_ones > n_zeroes) {
		
			diff <- n_ones - n_zeroes
			half_diff <- (diff / 2)
			swaps <- sample(ones, half_diff)
			ones <- ones[-which(ones %in% swaps)]
			zeroes <- c(zeroes, swaps)
		
		} else if (n_ones < n_zeroes) {
		
			diff <- n_zeroes - n_ones
			half_diff <- (diff / 2)
			swaps <- sample(zeroes, half_diff)
			zeroes <- zeroes[-which(zeroes %in% swaps)]
			ones <- c(ones, swaps)
		
		}
		response_predictors_no_loss$fold[ones] <- 1
		response_predictors_no_loss$fold[zeroes] <- 0

		# divide between calibration/evaluation sites
		calib <- rbind(
			response_predictors_loss[response_predictors_loss$fold == 0, ],
			response_predictors_no_loss[response_predictors_no_loss$fold == 0, ]
		)
		
		eval <- rbind(
			response_predictors_loss[response_predictors_loss$fold == 1, ],
			response_predictors_no_loss[response_predictors_no_loss$fold == 1, ]
		)

		n_calib <- nrow(calib)
		n_eval <- nrow(eval)

		n_calib_loss <- sum(calib$forest_loss)
		n_calib_no_loss <- sum(calib$forest_loss == 0)

		n_eval_loss <- sum(eval$forest_loss)
		n_eval_no_loss <- sum(calib$forest_loss == 0)

		### calibrate multivariate model
		################################

		form <- paste0('forest_loss ~ 1 + ', paste(predictors_selected, collapse = ' + '))
		form <- as.formula(form)
		
		if (any(predictors_selected == "forest_frag_class")) {
		
			calib$forest_frag_class <- factor(calib$forest_frag_class)
			eval$forest_frag_class <- factor(eval$forest_frag_class)

		}
		
		if (any(predictors_selected == "country")) {
		
			calib$country <- factor(calib$country)
			eval$country <- factor(eval$country)

		}

		model_multivar <- speedglm(form, data = calib, family = binomial())

		say('')
		print(summary(model_multivar))
		say('')

		# AUC
		predict_eval <- predict(model_multivar, eval, type = 'response')
		pred_eval <- prediction(predict_eval, labels = eval$forest_loss)
		auc_obs <- performance(pred_eval, 'auc')@y.values[[1]]
		
		# get coefficient values
		coefs <- coefficients(model_multivar)
		coef_names <- names(coefs)
		
		results <- rbind(
			results,
			data.frame(
				demesne = demesne,
				model_type = 'observed',
				k = k,
				n_calib = n_calib,
				n_eval = n_eval,
				n_calib_loss = n_calib_loss,
				n_calib_no_loss = n_calib_no_loss,
				n_eval_loss = n_eval_loss,
				n_eval_no_loss = n_eval_no_loss,
				coefficient = coef_names,
				coefficient_value = coefs,
				auc = auc_obs
			)
		)

		### univariate models and permutation test on multivariate model
		n <- nrow(eval)

		for (this_pred in predictors_selected) {

			### univariate model
			####################
			
			form <- paste0('forest_loss ~ 1 + ', this_pred)
			form <- as.formula(form)
			
			model_univar <- speedglm(form, data = calib, family = binomial())

			# AUC
			predict_eval <- predict(model_univar, eval, type = 'response')
			pred_eval <- prediction(predict_eval, labels = eval$forest_loss)
			auc_univar <- performance(pred_eval, 'auc')@y.values[[1]]
			
			# get coefficient values
			coef <- coefficients(model_univar)[this_pred]
			coef <- unname(coef)

			results <- rbind(
				results,
				data.frame(
					demesne = demesne,
					model_type = 'univariate',
					k = k,
					n_calib = n_calib,
					n_eval = n_eval,
					n_calib_loss = n_calib_loss,
					n_calib_no_loss = n_calib_no_loss,
					n_eval_loss = n_eval_loss,
					n_eval_no_loss = n_eval_no_loss,
					coefficient = this_pred,
					coefficient_value = coef,
					auc = auc_univar
				)
			)

			### leave-out permutation test
			##############################
			
			eval_copy <- eval
			auc_permute <- rep(NA_real_, n_permute)

			# iterate permutation test
			for (i in seq_len(n_permute)) {
			
				eval_copy[[this_pred]] <- sample(eval_copy[[this_pred]], n)
				predict_eval_permute <- predict(model_multivar, eval_copy, type = 'response')

				pred_eval <- prediction(predict_eval_permute, labels = eval$forest_loss)
				auc_permute[i] <- performance(pred_eval, 'auc')@y.values[[1]]
			
			}

			auc_permute <- mean(auc_permute)

			results <- rbind(
				results,
				data.frame(
					demesne = demesne,
					model_type = 'permuted',
					k = k,
					n_calib = n_calib,
					n_eval = n_eval,
					n_calib_loss = n_calib_loss,
					n_calib_no_loss = n_calib_no_loss,
					n_eval_loss = n_eval_loss,
					n_eval_no_loss = n_eval_no_loss,
					coefficient = this_pred,
					coefficient_value = NA,
					auc = auc_permute
				)
			)

		} # next predictor

		say('CROSS-VALIDATION FOLD ', k, ' WRITE PREDICTION RASTERS', level = 2)
		step <- paste0('Cross-validation: Write prediction rasters')
		
		# predict model raster
		# subset environmental rasters (makes prediction faster)
		start <- Sys.time()
		env_selected <- env[[predictors_selected]]
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_double_brackets', target = 'Predictor rasters', dtype = 'raster', start = start, stop = stop, n = length(env))
	
		# prediction raster
		start <- Sys.time()
		prediction_rast <- predict(env_selected, model_multivar, type = 'response')
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'predict()', target = 'Predictor raster', dtype = 'raster', start = start, stop = stop, n = length(env_selected))

		# save prediction raster
		fn <- paste0('./outputs/', tolower(demesne), '_fasterRaster_fold_', prefix(k, 2), '.tif')
		start <- Sys.time()
		writeRaster(prediction_rast, filename = fn, overwrite = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Predictor raster', dtype = 'raster', start = start, stop = stop)

	} # next fold
	
	write.csv(results, paste0('./outputs/', tolower(demesne), '_fasterRaster_results.csv'), row.names = FALSE)

	### POST-PROCESS PREDICTION RASTERS
	###################################
	
	say('POST-PROCESS PREDICTION RASTERS', level = 2)
	step <- 'Post-process prediction rasters'

	# load prediction rasters
	fns <- paste0('./outputs/', tolower(demesne), '_fasterRaster_fold_', prefix(1:n_folds, 2), '.tif')
	start <- Sys.time()
	predictions <- fast(fns)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictions))

	# average and range of rasters
	start <- Sys.time()
	prediction_mean <- mean(predictions)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mean()', target = 'Prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictions))

	start <- Sys.time()
	prediction_sd <- app(predictions, fun = 'sd', cores = cores)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'app()', target = 'Minimum & maximum prediction rasters', dtype = 'raster', start = start, stop = stop)

	# classify prediction rasters into 9 classes based on combination of probability of forest loss and uncertainty
	start <- Sys.time()
	quants_mean <- global(prediction_mean, fun = quantile, probs = threshold_quantiles, na.rm = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'global()', target = 'Mean prediction raster', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	quants_sd <- global(prediction_sd, fun = quantile, probs = threshold_quantiles, na.rm = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'quantile()', target = 'SD of predictions raster', dtype = 'raster', start = start, stop = stop)

	quants_mean <- unlist(quants_mean)
	quants_sd <- unlist(quants_sd)

	# stack mean and sd prediction
	start <- Sys.time()
	prediction_mean_sd <- c(prediction_mean, prediction_sd)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Mean and SD of predictions rasters', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(prediction_mean_sd) <- c('mu', 'sigma')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Mean and SD of predictions rasters', dtype = 'raster', start = start, stop = stop)

	# calculate classes based on low/medium/high probability of forest loss and uncertainty therein
	zone_fx <- function(mu, sigma, quants_mean, quants_sd) {
	
		mu_lower <- quants_mean[1]
		mu_upper <- quants_mean[2]

		sigma_lower <- quants_sd[1]
		sigma_upper <- quants_sd[2]

		out <- mu * NA
		out[mu <= mu_lower & sigma <= sigma_lower] <- 0 # low prob/low uncert
		out[mu > mu_lower & mu <= mu_upper & sigma <= sigma_lower] <- 1 # medium prob/low uncert
		out[mu > mu_upper & sigma <= sigma_lower] <- 2 # high prob/low uncert

		out[mu <= mu_lower & sigma > sigma_lower & sigma <= sigma_upper] <- 3 # low prob/medium uncert
		out[mu > mu_lower & mu <= mu_upper & sigma > sigma_lower & sigma <= sigma_upper] <- 4 # medium prob/medium uncert
		out[mu >= mu_upper & sigma > sigma_lower & sigma <= sigma_upper] <- 5 # high prob/medium uncert

		out[mu <= mu_lower & sigma > sigma_upper] <- 6 # low prob/high uncert
		out[mu > mu_lower & mu <= mu_upper & sigma > sigma_upper] <- 7 # medium prob/high uncert
		out[mu > mu_upper & sigma > sigma_upper] <- 8 # high prob/high uncert
		out

	}
	
	# NB using cores > 1 throws error: first error: "unused arguments (quants_mean = c(0.2665024 . . .))"
	start <- Sys.time()
	prediction_class <- lapp(prediction_mean_sd, fun = zone_fx, quants_mean = quants_mean, quants_sd = quants_sd, usenames = TRUE, cores = 1)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'lapp()', target = 'Mean and SD of predictions rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_sd))
	
	# name
	start <- Sys.time()
	names(prediction_class) <- 'prediction_class'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Predictions class raster', dtype = 'raster', start = start, stop = stop)

	# combine
	start <- Sys.time()
	prediction_mean_sd_class <- c(prediction_mean_sd, prediction_class)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Mean, SD, and class of prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_sd_class))

	# save prediction rasters (mean, sd, class)
	fn <- paste0('./outputs/', tolower(demesne), '_fasterRaster_prediction_mean_sd_class.tif')
	start <- Sys.time()
	writeRaster(prediction_mean_sd_class, fn, overwrite = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Mean, SD, and class of prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_sd_class))

	write.csv(timings, paste0('./outputs/', tolower(demesne), '_fasterRaster_timings.csv'), row.names = FALSE)

say('DONE!', deco = '^', level = 1)
say(date())
sink()
