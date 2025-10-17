### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, fasterRaster, and sf packages for R.
###
### source('C:/fasterRaster_benchmarking_distributed/fasterRaster_benchmarking/02_benchmarking_fasterRaster.r')
###
### CONTENTS ###
### settings ###
### benchmark ###

################
### settings ###
################

	rm(list=ls())

	# drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'
	drive <- 'C:/'

	.libPaths(paste0(drive, '/fasterRaster_benchmarking_distributed/libraries'))

	# setwd(paste0(drive, '/Research/fasterRaster - Streamlined GIS in R through GRASS'))
	setwd(paste0(drive, '/fasterRaster_benchmarking_distributed'))

	source('./fasterRaster_benchmarking/00_constants.r')

	memory <- 64 * 1024
	# memory <- 128 * 1024 # if not comparing to terra
	cores <- 4

	grass_dir <- 'C:/Program Files/GRASS GIS 8.4/'
	faster(grassDir = grass_dir, memory = memory, cores = cores, useDataTable = TRUE, verbose = TRUE)

	threads <- TRUE
	memmax <- 64
	# memmax <- 128 # if not comparing to terra

	terraOptions(
		memmax = memmax,
		progress = 0
	)
	
	x11() # start diplay bc later will plot slope for manual inspection of anomolies

say('#################')
say('### benchmark ###')
say('#################')

	### options
	###########

	# demesne <- 'Large'
	# demesne <- 'Medium'
	demesne <- 'Small'

	### start
	#########

	if (demesne == 'Small') {

		basins_primary <- 'Mekong'
		basins_secondary <- 'Nam Loi'
		country_names <- c('China', 'Myanmar')
		n_folds <- n_folds_small
		glad_lab_tile_names <- c('30N_090E', '30N_100E')
		pop_tile_names <- 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C28'
		extent <- c(73, 134, 9, 54)
		aggregate_factor <- small_agg_factor
		cross_valid_n <- cross_valid_n_small

	} else if (demesne == 'Medium') {

		basins_primary <- 'Salween'
		basins_secondary <- NA
		country_names <- c('China', 'Myanmar', 'Thailand')
		n_folds <- n_folds_medium
		glad_lab_tile_names <- c('30N_090E', '40N_090E', '20N_090E', '30N_100E')
		pop_tile_names <- c(
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C27',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C28',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C28',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R8_C28'
		)
		extent <- c(73, 134, 8, 54)
		aggregate_factor <- medium_agg_factor
		cross_valid_n <- cross_valid_n_medium

	} else if (demesne == 'Large') {

		basins_primary <- c('Mekong', 'Salween', 'Irrawaddy', 'Chao Phraya', 'Sittang')
		basins_secondary <- NA
		country_names <- c('China', 'India', 'Myanmar', 'Thailand', 'Laos', 'Cambodia', 'Vietnam')
		n_folds <- n_folds_large
		glad_lab_tile_names <- c('30N_090E', '20N_100E', '20N_090E', '10N_100E', '40N_090E', '30N_100E')
		pop_tile_names <- c(
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R5_C27', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C27',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C27',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C29',
			'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R8_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R8_C29'
		)
		extent <- c(73, 135, 8, 54)
		aggregate_factor <- large_agg_factor
		cross_valid_n <- cross_valid_n_large

	}

	# output_dir <- paste0(substr(drive, 1, 2), '/!scratch/terra_outputs_', tolower(demesne), '/')
	# output_dir <- paste0('./fasterRaster_outputs_', tolower(demesne), '/', ifelse(aggregate_factor == 1, '_dist_agg_factor_1', ''))
	output_dir <- paste0('./fasterRaster_outputs_', tolower(demesne), '_agg_factor_', aggregate_factor, '_attempt/')
	dirCreate(output_dir)

	sink(paste0(output_dir, tolower(demesne), '_fasterRaster_benchmark.txt'), split = TRUE)
	say('BENCHMARK FASTERRASTER')
	say(date(), post = 1)

	say('SESSION INFO', level = 2)
	print(sessionInfo())

	say('FASTERRASTER OPTIONS', level = 2)
	print(faster())

	say('TERRA OPTIONS', level = 2)
	print(terraOptions())

	say('DEMESNE', level = 2)
	say('Demesne: ....................................................... ', demesne)
	say('Primary Basins: ................................................ ', paste(basins_primary, collapse = ', '))
	say('Secondary Basins: .............................................. ', paste(basins_secondary, collapse = ', '))
	say('Country Names: ................................................. ', paste(country_names, collapse = ', '), post = 1)

	say('SETTINGS', level = 2)
	say('Drive with datasets and libraries: ............................... ', drive)
	say('Number of cores used for multi-core functions: ................... ', cores)
	say('Use multiple threads for `terra` functions that allow it: ........ ', threads)
	say('Maximum memory allowed for fasterRaster (GB): .................... ', memory / 1024)
	say('Maximum memory allowed for terra (GB): ........................... ', memmax)
	say('Study region buffer size (m): .................................... ', study_region_buffer_size_m)
	say('Factor by which to aggregate raster for distance calculations: ... ', aggregate_factor)
	say('Number of folds: ................................................. ', n_folds)
	say('Number of sites per fold: ........................................ ', cross_valid_n)
	say('Inflation of number of points for variable selection: ............ ', inflation_for_variable_selection)
	say('Quantiles used to divide final prediction raster into zones: ..... ', paste(threshold_quantiles, collapse = ', '))
	say('Number of iterations in permutation importance test: ............. ', n_permute)
	say('')

	all_ok <- function() {
		x <- readline('Are the settings above OK? (y/n) ')
		if (x != 'y') stop('All your base are belong to us.')
	}
	
	all_ok()

	all_ok <- function() {
		x <- readline('Have you already done a `FINAL` run using these parameters (y/n)? ')
		if (x != 'n') stop('All your base are belong to us.')
	}
	
	all_ok()

	# stores runtime information
	timings <- data.table()

	### LOAD RAW RASTERS
	####################

	say('LOAD RAW RASTERS', level = 2)
	step <- 'Load raw rasters'

	# # forest cover for 2000 and 2020
	# forest_files_2000 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2020 Forest Height & Extent/2000/', glad_lab_tile_names, '.tif')
	# forest_files_2020 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2020 Forest Height & Extent/2020/', glad_lab_tile_names, '.tif')
	
	# forest cover for 2000 and 2020
	forest_files_2000 <- paste0('./Potapov et al 2020 Forest Height & Extent/2000/', glad_lab_tile_names, '.tif')
	forest_files_2020 <- paste0('./Potapov et al 2020 Forest Height & Extent/2020/', glad_lab_tile_names, '.tif')
	
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
	
	# # LULC
	# lulc_files_2000 <- paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2022 LULC/2000/', glad_lab_tile_names, '.tif')
	# lulc_levels <- read.csv(paste0(drive, '/Research Data/GLAD Lab/Potapov et al 2022 LULC/legend.csv'))

	lulc_files_2000 <- paste0('./Potapov et al 2022 LULC/2000/', glad_lab_tile_names, '.tif')
	lulc_levels <- read.csv(paste0('./Potapov et al 2022 LULC/legend.csv'))

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
	# start <- Sys.time()
	# population <- fast(paste0(drive, '/Research Data/Population - Gridded Population of the World (GPW) v4/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_30_sec.tif'))
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'fast()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	population <- list()
	start <- Sys.time()
	for (i in seq_along(pop_tile_names)) {

		pop_file <- paste0('./Human Population - GHSL/R2023_Mollweide/', pop_tile_names[i], '/', pop_tile_names[i], '.tif')
		population[[i]] <- fast(pop_file)

	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Population density', dtype = 'raster', start = start, stop = stop, n = length(pop_tile_names))


	# # soil nitrogen
	# start <- Sys.time()
	# nitrogen <- fast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean.tif'))
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'fast()', target = 'Soil nitrogen', dtype = 'raster', start = start, stop = stop)

	# # soil organic carbon
	# start <- Sys.time()
	# soc <- fast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean.tif'))
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'fast()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	### LOAD RAW VECTORS I
	######################

	say('LOAD RAW VECTORS I', level = 2)
	step <- 'Load raw vectors'

	# drainage basins
	start <- Sys.time()

	# se_asia_basins <- fast(paste0(drive, '/Research Data/Watershed Basins in Southeast Asia (FAO)/hydrobasins_asia.gpkg'), extent = extent)
	se_asia_basins <- fast('./Watershed Basins in Southeast Asia (FAO)/hydrobasins_asia.gpkg', extent = extent)

	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Watershed basins', dtype = 'vector', start = start, stop = stop)

	### DEFINE STUDY REGION
	#######################

	say('DEFINE STUDY REGION', level = 2)
	step <- 'Define study region'
	
	### select study region
	start <- Sys.time()
	if (demesne == 'Small') {
		focal_basins <- se_asia_basins[se_asia_basins$MAJ_NAME %in% basins_primary & se_asia_basins$SUB_NAME %in% basins_secondary]
	} else {
		focal_basins <- se_asia_basins[se_asia_basins$MAJ_NAME %in% basins_primary]
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Basins vector', dtype = 'vector', start = start, stop = stop)

	# aggregate polygons of basins to create study region
	# we don't need to do this for the small region because it's just one basin, but we do it anyway to assess comparability
	start <- Sys.time()
	study_region <- aggregate(focal_basins)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'aggregate()', target = 'Mekong basin vector', dtype = 'vector', start = start, stop = stop)

	# remove holes
	if (grassInfo('versionNumber') >= 8.4) {
		
		start <- Sys.time()
		study_region <- fillHoles(study_region)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'fillHoles()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	}

	# remove 'islands'
	start <- Sys.time()
	study_region <- disagg(study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'disagg()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	start <- Sys.time()
	study_region_areas_m2 <- expanse(study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'expanse()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	start <- Sys.time()
	study_region <- study_region[which.max(study_region_areas_m2)]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	# add buffer to basins to obviate edge effects
	# For unprojected CRSs, GRASS buffer sizes are in degrees. So we estimate the desire buffer size using the mean longitude of the study region.
	mean_lat_deg <- mean(ext(focal_basins, TRUE)[3:4])
	mean_lat_rad <- (mean_lat_deg / 90) * (pi / 2)
	study_region_buffer_size_deg <- cos(mean_lat_rad) * study_region_buffer_size_m / 111120

	start <- Sys.time()
	study_region_buffered <- buffer(study_region, dissolve = TRUE, width = study_region_buffer_size_deg) # buffer in degrees
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'buffer()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	# save study region vector
	start <- Sys.time()
	writeVector(study_region, paste0(output_dir, tolower(demesne), '_fasterRaster_study_region.gpkg'), overwrite = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeVector()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	### LOAD RAW VECTORS II
	#######################

	say('LOAD RAW VECTORS II', level = 2)
	step <- 'Load raw vectors'

	# rivers
	start <- Sys.time()
	# rivers <- fast(paste0(drive, '/Research Data/Rivers - FAO/rivers_asia_37331.shp'), extent = study_region_buffered)
	rivers <- fast('./Rivers - FAO/rivers_asia_37331.shp', extent = study_region_buffered)
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
		
		# fn <- paste0(drive, '/Research Data/Open Street Maps (Geofabrik)/', country, '-latest-free.shp/gis_osm_roads_free_1.shp')
		fn <- paste0('./Open Street Maps (Geofabrik)/', country, '-latest-free.shp/gis_osm_roads_free_1.shp')
		
		roads[[i]] <- fast(fn, dropTable = TRUE, extent = study_region_buffered)
	
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'Roads', dtype = 'vector', start = start, stop = stop, n = length(roads))

	# protected areas
	pas <- list()
	start <- Sys.time()
	for (i in 1:3) {

		# fn <- paste0(drive, '/Research Data/World Database of Protected Areas/WDPA_May2024_Public_shp_', i - 1, '/WDPA_May2024_Public_shp-polygons.shp')
		fn <- paste0('./World Database of Protected Areas/WDPA_May2024_Public_shp_', i - 1, '/WDPA_May2024_Public_shp-polygons.shp')

		pas[[length(pas) + 1]] <- fast(fn, extent = study_region_buffered, resolve = 'aggregate')
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'WDPA', dtype = 'vector', start = start, stop = stop, n = 3)

	# countries of SE Asia
	start <- Sys.time()
	# countries <- fast(paste0(drive, '/Research Data/GADM/Version 4.1/gadm_410.gpkg'), extent = study_region_buffered)
	countries <- fast('./GADM/Version 4.1/gadm_410.gpkg', extent = study_region_buffered)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'fast()', target = 'SE Asian countries', dtype = 'vector', start = start, stop = stop)

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
	
	# rasterize study region to use as mask
	start <- Sys.time()
	study_region_rast <- rasterize(study_region_buffered, forest_2000)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rasterize()', target = 'Study region', dtype = 'vector/raster', start = start, stop = stop)

	# trim study region
	start <- Sys.time()
	study_region_rast <- trim(study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'trim()', target = 'Study region raster', dtype = 'raster', start = start, stop = stop)

	# crop
	start <- Sys.time()
	forest_2000 <- crop(forest_2000, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = '2000 forest cover', dtype = 'raster', start = start, stop = stop)
	
	start <- Sys.time()
	forest_2020 <- crop(forest_2020, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = '2020 forest cover', dtype = 'raster', start = start, stop = stop)
	
	# aggregated study region
	# aggregate raster to make processing faster
	start <- Sys.time()
	study_region_rast_agg <- aggregate(study_region_rast, fact = aggregate_factor)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'aggregate()', target = 'Study region raster', dtype = 'raster', start = start, stop = stop)
	
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
	forest_2000_no_forest_masked <- forest_2000

	start <- Sys.time()
	forest_2000_no_forest_masked[forest_2000_no_forest_masked == 0] <- NA_integer_
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'replace_single_square_bracket', target = 'Forest 2000 raster', dtype = 'raster', start = start, stop = stop)
	
	start <- Sys.time()
	forest_loss <- forest_loss * forest_2000_no_forest_masked
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Masked forest 2000 & loss raster', dtype = 'raster', start = start, stop = stop)
	
	# name
	start <- Sys.time()
	names(forest_loss) <- 'forest_loss'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)

	### SLOPE
	#########
	
	say('SLOPE', level = 2)
	step <- 'Prepare slope predictor'

# save and reload to obviate random horizontal band problem
dirCreate(paste0(output_dir, '/!scratch'))
writeRaster(elev_scale_7_m, paste0(output_dir, '/!scratch/', tolower(demesne), '_elev_scale_7_m.tif'), overwrite = TRUE, bigTiff = TRUE)
elev_scale_7_m <- fast(paste0(output_dir, '/!scratch/', tolower(demesne), '_elev_scale_7_m.tif'))

	# slope
	start <- Sys.time()
	slope_scale_7 <- terrain(elev_scale_7_m, 'slope')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'terrain()', target = 'Elevation at scale 7', dtype = 'raster', start = start, stop = stop)

# save and reload to obviate random horizontal band problem
writeRaster(slope_scale_7, paste0(output_dir, '/!scratch/', tolower(demesne), '_slope_scale_7_STEP_1.tif'), overwrite = TRUE, bigTiff = TRUE)
slope_scale_7 <- fast(paste0(output_dir, '/!scratch/', tolower(demesne), '_slope_scale_7_STEP_1.tif'))

	start <- Sys.time()
	slope_scale_11 <- terrain(elev_scale_11_m, 'slope')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'terrain()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	# resample
	start <- Sys.time()
	slope_scale_7 <- resample(slope_scale_7, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Slope at scale 7', dtype = 'raster', start = start, stop = stop)

# save and reload to obviate random horizontal band problem... doesn't work!
writeRaster(slope_scale_7, paste0(output_dir, '/!scratch/', tolower(demesne), '_slope_scale_7_STEP_2.tif'), overwrite = TRUE, bigTiff = TRUE)
slope_scale_7 <- fast(paste0(output_dir, '/!scratch/', tolower(demesne), '_slope_scale_7_STEP_2.tif'))

	start <- Sys.time()
	slope_scale_11 <- resample(slope_scale_11, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Slope at scale 11', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(slope_scale_7) <- 'slope_scale_7'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Slope at scale 7', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	names(slope_scale_11) <- 'slope_scale_11'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Slope at scale 11', dtype = 'raster', start = start, stop = stop)

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
	elev_scale_11_m <- resample(elev_scale_11_m, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	# remove impossibly negative value... artifact!?!
	# if terra does not have this issue, will have to retain this for overall performance but not include in fx-vs-fx comparison
	start <- Sys.time()
	impossible <- minmax(elev_scale_11_m)['min', ] < 1
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'minmax()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)

	if (impossible) {
	
		start <- Sys.time()
		elev_scale_11_m[elev_scale_11_m < 0] <- NA
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'replace_single_square_bracket', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)
	
	}

	# names
	start <- Sys.time()
	names(elev_scale_11_m) <- 'elev_scale_11_m'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Elevation at scale 11', dtype = 'raster', start = start, stop = stop)
	
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
	dist_to_rivers_km <- resample(dist_to_rivers_agg_km, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = '2000 forest cover (coarse and fine scale)', dtype = 'raster', start = start, stop = stop)

	# log
	start <- Sys.time()
	dist_to_rivers_km <- dist_to_rivers_km + 1
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = '2000 forest cover (coarse and fine scale)', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	log10_dist_to_rivers_km <- log10(dist_to_rivers_km)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'log10()', target = '2000 forest cover (coarse and fine scale)', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(log10_dist_to_rivers_km) <- 'log10_dist_to_rivers_km'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Distance to rivers', dtype = 'raster', start = start, stop = stop)
	
	### COUNTRIES
	#############
	
	say('COUNTRIES', level = 2)
	step <- 'Prepare country predictor'

	# rasterize
	start <- Sys.time()
	countries_rast <- rasterize(countries, study_region_rast, field = 'COUNTRY')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rasterize()', target = 'Countries vector', dtype = 'vector/raster', start = start, stop = stop)
	
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
	timings <- remember(timings, step = step, fx = 'merge()', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = length(glad_lab_tile_names))

	# crop
	start <- Sys.time()
	lulc <- crop(lulc, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'LULC', dtype = 'raster', start = start, stop = stop)

	# mask short vegetation
	start <- Sys.time()
	short_veg <- (lulc >= 19 & lulc <= 24) | (lulc >= 102 & lulc <= 124)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'comparison', target = 'LULC', dtype = 'raster', start = start, stop = stop, n = 4)

	# density of short vegetation @ scale 2
	start <- Sys.time()
	short_veg_density_33cells <- focal(short_veg, w = 33, fun = 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'focal()', target = 'Short vegetation', dtype = 'raster', start = start, stop = stop)
	
	# name
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
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Roads', dtype = 'raster', start = start, stop = stop, n = length(roads))
	
	# combine roads vectors
	start <- Sys.time()
	roads <- do.call(rbind, roads)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rbind()', target = 'Roads', dtype = 'raster', start = start, stop = stop, n = length(country_names))
	
	# distance to nearest road
	start <- Sys.time()
	dist_to_roads_agg_km <- distance(study_region_rast_agg, roads, unit = 'km')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'distance()', target = 'Aggregated study region raster & roads', dtype = 'raster/vector', start = start, stop = stop)
	
	# resample
	start <- Sys.time()
	dist_to_roads_km <- resample(dist_to_roads_agg_km, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'resample()', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)

	# log
	start <- Sys.time()
	dist_to_roads_km <- dist_to_roads_km + 1
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	log10_dist_to_roads_km <- log10(dist_to_roads_km)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'log10()', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(log10_dist_to_roads_km) <- 'log10_dist_to_roads_km'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Distance to roads', dtype = 'raster', start = start, stop = stop)
	
	# ### SOIL
	# ########
	
	# say('SOIL', level = 2)
	# step <- 'Prepare soil predictors'

	# # project study region to soil CRS
	# start <- Sys.time()
	# study_region_igh <- project(study_region_buffered, nitrogen)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'project()', target = 'Study region', dtype = 'vector', start = start, stop = stop)

	# # crop soil rasters to study region
	# start <- Sys.time()
	# nitrogen <- crop(nitrogen, study_region_igh)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'crop()', target = 'Nitrogen', dtype = 'raster', start = start, stop = stop)

	# start <- Sys.time()
	# soc <- crop(soc, study_region_igh)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'crop()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	# # project soil rasters to study region CRS and forest resolution
	# start <- Sys.time()
	# nitrogen <- project(nitrogen, study_region_rast, method = 'bilinear')
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'project()', target = 'Nitrogen', dtype = 'raster', start = start, stop = stop)

	# start <- Sys.time()
	# soc <- project(soc, study_region_rast, method = 'bilinear')
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'project()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	# # crop soil rasters to study region
	# start <- Sys.time()
	# nitrogen <- crop(nitrogen, study_region_rast)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'crop()', target = 'Nitrogen', dtype = 'raster', start = start, stop = stop)

	# start <- Sys.time()
	# soc <- crop(soc, study_region_rast)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'crop()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	# # names
	# start <- Sys.time()
	# names(nitrogen) <- 'nitrogen'
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'names()', target = 'Nitrogen', dtype = 'raster', start = start, stop = stop)

	# start <- Sys.time()
	# names(soc) <- 'soc'
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'names()', target = 'Soil organic carbon', dtype = 'raster', start = start, stop = stop)

	### HUMAN POPULATION
	####################

	say('HUMAN POPULATION', level = 2)
	step <- 'Prepare human population density predictor'

	if (demesne %in% c('Large', 'Medium')) {
		
		# merge
		start <- Sys.time()
		population <- do.call(merge, population)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'merge()', target = 'Population density', dtype = 'raster', start = start, stop = stop, n = length(pop_tile_names))

	} else {
		population <- population[[1]]
	}

	# project
	start <- Sys.time()
	population <- project(population, study_region_rast, method = 'bilinear')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'project()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	# crop
	start <- Sys.time()
	population <- crop(population, study_region_rast)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Population density', dtype = 'raster', start = start, stop = stop)
	
	# focal
	start <- Sys.time()
	population <- focal(population, w = 33, fun = 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	# transform
	start <- Sys.time()
	population <- population + 1
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	population <- log10(population)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'log10()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(population) <- 'log10_population'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Population density', dtype = 'raster', start = start, stop = stop)

	### PROTECTED AREAS
	###################

	say('PROTECTED AREAS', level = 2)
	step <- 'Prepare protected areas predictor'
	
	# crop
	pas_study_region <- list()
	start <- Sys.time()
	for (i in seq_along(pas)) {
		pas_study_region[[length(pas_study_region) + 1]] <- crop(pas[[i]], study_region_rast, fail = FALSE)
	}
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Protected areas', dtype = 'vector', start = start, stop = stop, n = length(pas))

	# combine
	if (length(pas_study_region) > 1) {
		
		start <- Sys.time()
		pas <- do.call(rbind, pas_study_region)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'rbind()', target = 'Protected areas', dtype = 'vector', start = start, stop = stop, n = length(pas))

	} else {
		pas <- pas_study_region[[1]]
	}

	# rasterize
	start <- Sys.time()
	pas_rast <- rasterize(pas, study_region_rast, background = 0)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'rasterize()', target = 'Protected areas', dtype = 'vector/raster', start = start, stop = stop)

	# add levels
	start <- Sys.time()
	levels(pas_rast) <- data.table(value = 0L:1L, label = c('unprotected', 'protected'))
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'levels()', target = 'Protected areas', dtype = 'raster', start = start, stop = stop)

	# names
	start <- Sys.time()
	names(pas_rast) <- 'PA'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Protected areas', dtype = 'raster', start = start, stop = stop)

	### COLLATE PREDICTOR RASTERS
	#############################

	say('COLLATE PREDICTOR RASTERS', level = 2)
	step <- 'Collate predictor rasters'

	### stack continuous predictors
	start <- Sys.time()
	env <- c(
		forest_loss,
		forest_frag_class,
		countries_rast,
		pas_rast,
		forest_density_33cells,
		elev_scale_11_m,
		slope_scale_7,
		slope_scale_11,
		log10_dist_to_rivers_km,
		log10_dist_to_roads_km,
		short_veg_density_33cells,
		# nitrogen,
		# soc,
		population
	)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	# plot for checking
	plot(env$slope_scale_7, main = 'slope_scale_7 from fR', simplify = FALSE)

	# remove India if analyzing large region because too little protected area there to be statistically stable
	if (demesne == 'Large') {
	
		# levels
		start <- Sys.time()
		country_levels <- levels(countries_rast)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'levels()', target = 'Countries raster', dtype = 'raster', start = start, stop = stop)

		india <- country_levels[[1]]$value[which(country_levels[[1]]$COUNTRY == 'India')]

		# names
		start <- Sys.time()
		env_names <- names(env)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'names()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

		# mask
		start <- Sys.time()
		env <- mask(env, mask = countries_rast, maskvalues = india, inverse = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'mask()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

		# names
		start <- Sys.time()
		names(env) <- env_names
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'names()<-', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	}

	# crop and mask to non-buffered study region
	start <- Sys.time()
	env <- crop(env, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'crop()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	start <- Sys.time()
	env <- mask(env, study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Response and continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	### make mask across rasters where any non-NA cell is 1 and any NA cell is NA
	start <- Sys.time()
	na_nonna_mask <- sum(env)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'sum()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	start <- Sys.time()
	env <- mask(env, na_nonna_mask)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mask()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	start <- Sys.time()
	env <- trim(env)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'trim()', target = 'Continuous & non-continuous rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	### save response and predictor rasters
	bigTiff <- demesne %in% c('Large', 'Medium')
	start <- Sys.time()
	writeRaster(env, paste0(output_dir, tolower(demesne), '_fasterRaster_response_predictors.tif'), overwrite = TRUE, bigTiff = bigTiff, byLayer = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	# retain forest loss raster for selecting model calibration/evaluation sites
	start <- Sys.time()
	forest_loss_focal_basin <- env[['forest_loss']]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_brackets', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)
	
	### scale continuous predictors
	# names_continuous_predictors <- c('forest_density_33cells', 'elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'short_veg_density_33cells', 'nitrogen', 'soc', 'log10_population')
	names_continuous_predictors <- c('forest_density_33cells', 'elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'short_veg_density_33cells', 'log10_population')

	start <- Sys.time()	
	continuous_predictors <- env[[names_continuous_predictors]]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_brackets', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(continuous_predictors))
		
	start <- Sys.time()
	continuous_predictors_scaled <- scalepop(continuous_predictors)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'scale()', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(continuous_predictors))

	start <- Sys.time()
	env <- env[[!(names(env) %in% names_continuous_predictors)]]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_brackets', target = 'Continuous predictor rasters', dtype = 'raster', start = start, stop = stop, n = length(names_continuous_predictors))

	start <- Sys.time()
	env <- c(env, continuous_predictors_scaled)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	### save scaled response and predictor rasters
	bigTiff <- demesne %in% c('Large', 'Medium')
	start <- Sys.time()
	writeRaster(env, paste0(output_dir, tolower(demesne), '_scaled_fasterRaster_response_predictors.tif'), overwrite = TRUE, bigTiff = bigTiff, byLayer = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))
	
	### CLEAN-UP DISK SPACE
	#######################
	
	say('REMOVE TEMPORARY FILES', level = 2)
	step <- 'Remove temporary files'

	keep <- list(continuous_predictors_scaled, forest_loss_focal_basin, env, study_region)
	start <- Sys.time()
	mow('*', ask = FALSE, keep = keep)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mow()', target = 'Temporary files', dtype = 'raster/vector', start = start, stop = stop)

	### CROSS-VALIDATION SETUP
	##########################
	
	say('CROSS-VALIDATION SETUP', level = 2)
	step <- 'Cross-validation setup'

	# calculate number of sites for calibration and evaluation
	# calculate sample size... selecting more than needed bc sites can be placed in NA cells

	start <- Sys.time()
	not_na <- not.na(forest_loss_focal_basin)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'not.na()', target = 'Forest loss raster', dtype = 'raster', start = start, stop = stop)

	# NB in fasterRaster, we could use nonnacell() to get this, but the underlying algorithm is the same
	start <- Sys.time()
	non_na_cells <- global(not_na, 'sum')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'global()', target = 'Forest loss/no loss raster', dtype = 'raster', start = start, stop = stop)

	non_na_cells <- non_na_cells$sum

	start <- Sys.time()
	cell_res <- res(forest_loss_focal_basin)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'res()', target = 'Forest loss/no loss raster', dtype = 'raster', start = start, stop = stop)

	non_na_cells_area_m2 <- prod(100000 * cell_res) * non_na_cells

	start <- Sys.time()
	study_region_area_m2 <- expanse(study_region)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'expanse()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	# number of random points to draw
	# inflating by 1.2 to ensure we have enough points on non-NA cells
	# also inflating by proportion of NA cells bc we will need to discard points that fall in these cells

	# fold sample size
	inflation <- 1.2 * study_region_area_m2 / non_na_cells_area_m2
	cross_valid_inflated_n <- ceiling(inflation * cross_valid_n)

	# variable selection sample size
	# inflated by user-specifed amount bc we want to ensure quasi-stable correlations
	# also inflating by 1.2 and ratio of non-NA to NA cells in raster
	variable_selection_inflated_n <- ceiling(inflation * inflation_for_variable_selection * cross_valid_n)
	variable_selection_n <- ceiling(inflation_for_variable_selection * cross_valid_n)

	say('non_na_cells ................... ', non_na_cells)
	say('initially want ', variable_selection_inflated_n ,' points for variable selection')
	say('eventually want ', variable_selection_n ,' points for variable selection')
	say('for each fold, initially selecting ', cross_valid_inflated_n ,' points, which will be subset to ', cross_valid_n, ' after removing NAs')

	### VARIABLE SELECTION
	######################

	say('VARIABLE SELECTION', level = 2)
	step <- 'Variable selection'

	start <- Sys.time()
	variable_selection_points <- spatSample(study_region, variable_selection_inflated_n, values = FALSE, as.points = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)

	start <- Sys.time()
	variable_selection_values <- extract(env, variable_selection_points)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'extract()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	# variable_selection_values <- as.data.table(variable_selection_values) # not needed for fasterRaster bc using data.table outputs
	variable_selection_values <- variable_selection_values[complete.cases(variable_selection_values), ]
	variable_selection_values <- variable_selection_values[sample(nrow(variable_selection_values), variable_selection_n), ]
	variable_selection_values <- variable_selection_values[ , c('forest_loss', 'forest_frag_class', 'country', 'PA') := NULL]

	### select variables with low collinearity
	cors <- vifcor(variable_selection_values, method = 'spearman', size = variable_selection_n, th = 0.7)

	say('Collinearity analysis and variable selection:', pre = 1)
	print(cors)
	say('')
	
	continuous_predictors_selected <- cors@results$Variables

	rm(cors)

	vars_selected <- c('forest_loss', continuous_predictors_selected, 'forest_frag_class', 'country', 'PA')
	predictors_selected <- c(continuous_predictors_selected, 'forest_frag_class', 'country', 'PA')

	start <- Sys.time()
	env <- env[[vars_selected]]
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'subset_double_brackets', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	### CROSS VALIDATION
	####################

	# NB We do not select random points from the full predictor stack first. Rather, we select random points from the forest loss/no loss raster, then use it to remove points that are outside the study region (NA cells) and to select the desired number of cells of each class (forest loss/persistence). *Then* we extract from the full response/predictor stack. This is much faster than selecting points and extracting values using the full stack initially.

	# store results from modeling in this data frame
	results <- data.table()
	prediction_rasts <- list()

	### CLEAN-UP DISK SPACE
	#######################
	
	say('REMOVE TEMPORARY FILES', level = 2)
	step <- 'Remove temporary files'

	keep <- list(study_region, env, forest_loss_focal_basin)

	start <- Sys.time()
	mow('*', ask = FALSE, keep = keep)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mow()', target = 'Temporary files', dtype = 'raster/vector', start = start, stop = stop)

	write.csv(data.frame(predictor = predictors_selected), paste0(output_dir, tolower(demesne), '_fasterRaster_predictors_selected.csv'), row.names = FALSE)
	fwrite(timings, paste0(output_dir, tolower(demesne), '_fasterRaster_timings_preliminary.csv'), row.names = FALSE)
	say('DONE with PRE-MODELING WORKFLOW!', deco = '^', level = 1)
	say(date())
# sink()
# print(NON)

	k <- 1
	while (k <= n_folds) {

		say('CROSS-VALIDATION FOLD ', k, level = 2)
		step <- paste0('Cross-validation: Select and extract from calibration and evaluation points')

		### randomly select points
		start <- Sys.time()
		sample_sites <- spatSample(study_region, cross_valid_inflated_n, as.points = TRUE, values = FALSE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop, k = k)

		start <- Sys.time()
		values_at_points <- extract(forest_loss_focal_basin, sample_sites)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'extract()', target = 'Forest loss raster masked to focal basin', dtype = 'raster', start = start, stop = stop, k = k)

		# remove points in NA cells
		completes <- which(complete.cases(values_at_points))
		start <- Sys.time()
		sample_sites <- sample_sites[completes]
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop, k = k)

		# keep desired number of point
		start <- Sys.time()
		n_sample_sites <- ngeom(sample_sites)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'ngeom()/nrow()', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop, k = k)

		keeps <- sample(n_sample_sites, cross_valid_n)

		start <- Sys.time()
		sample_sites <- sample_sites[keeps]
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop, k = k)

		### extract from full response/predictor stack
		start <- Sys.time()
		response_predictors <- extract(env, sample_sites, cats = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'extract()', target = 'Predictor and response rasters & calibration/evaluation sites', dtype = 'raster/vector', start = start, stop = stop, n = length(env), k = k)

		# catch cases with NAs
		nas <- !complete.cases(response_predictors)
		if (any(nas)) stop('NAs in `response_predictors`.')

		# divide into calibration/evaluation
		# Split loss/no loss into their own data frames, then assign random folds, then equilibrate so we have an equal number of cases in each fold in each loss/no loss group.
		response_predictors$fold <- NA_integer_		
		response_predictors_loss <- response_predictors[response_predictors$forest_loss == 1, ]
		response_predictors_no_loss <- response_predictors[response_predictors$forest_loss == 0, ]
		
		# apportion folds equally across LOSS sites
		n <- nrow(response_predictors_loss)
		fold <- sample(0L:1L, n, replace = TRUE)
		ones <- which(fold == 1)
		zeros <- which(fold == 0)
		
		n_ones <- length(ones)
		n_zeroes <- length(zeros)

		if (n_ones > n_zeroes) {
		
			diff <- n_ones - n_zeroes
			half_diff <- (diff / 2)
			swaps <- sample(ones, half_diff)
			ones <- ones[-which(ones %in% swaps)]
			zeros <- c(zeros, swaps)
		
		} else if (n_ones < n_zeroes) {
		
			diff <- n_zeroes - n_ones
			half_diff <- (diff / 2)
			swaps <- sample(zeros, half_diff)
			zeros <- zeros[-which(zeros %in% swaps)]
			ones <- c(ones, swaps)
		
		}
		response_predictors_loss$fold[ones] <- 1
		response_predictors_loss$fold[zeros] <- 0

		# apportion folds equally across NO LOSS sites
		n <- nrow(response_predictors_no_loss)
		fold <- sample(0L:1L, n, replace = TRUE)
		ones <- which(fold == 1)
		zeros <- which(fold == 0)
		
		n_ones <- length(ones)
		n_zeroes <- length(zeros)

		if (n_ones > n_zeroes) {
		
			diff <- n_ones - n_zeroes
			half_diff <- (diff / 2)
			swaps <- sample(ones, half_diff)
			ones <- ones[-which(ones %in% swaps)]
			zeros <- c(zeros, swaps)
		
		} else if (n_ones < n_zeroes) {
		
			diff <- n_zeroes - n_ones
			half_diff <- (diff / 2)
			swaps <- sample(zeros, half_diff)
			zeros <- zeros[-which(zeros %in% swaps)]
			ones <- c(ones, swaps)
		
		}
		response_predictors_no_loss$fold[ones] <- 1
		response_predictors_no_loss$fold[zeros] <- 0

		# divide between calibration/evaluation sites
		calib <- rbind(
			response_predictors_loss[fold == 0, ],
			response_predictors_no_loss[fold == 0, ]
		)
		
		eval <- rbind(
			response_predictors_loss[fold == 1, ],
			response_predictors_no_loss[fold == 1, ]
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
		if (demesne != 'Small') form <- paste0(form, ' + country:PA') # add interaction between country and protected area

		if (any(predictors_selected == 'log10_population')) form <- paste0(form, ' + I(log10_population^2)')
		if (any(predictors_selected == 'log10_dist_to_roads_km')) form <- paste0(form, ' + I(log10_dist_to_roads_km^2)')
		if (any(predictors_selected == 'log10_dist_to_rivers_km')) form <- paste0(form, ' + I(log10_dist_to_rivers_km^2)')
		if (any(predictors_selected == 'short_veg_density_33cells')) form <- paste0(form, ' + I(short_veg_density_33cells^2)')
		if (any(predictors_selected == 'forest_density_33cells')) form <- paste0(form, ' + I(forest_density_33cells^2)')
		
		form <- as.formula(form)

		if (any(predictors_selected == 'forest_frag_class')) {
		
			calib$forest_frag_class <- factor(calib$forest_frag_class)
			eval$forest_frag_class <- factor(eval$forest_frag_class)

		}
		
		if (any(predictors_selected == 'country')) {
		
			calib$country <- factor(calib$country)
			eval$country <- factor(eval$country)

		}

		# model_multivar <- speedglm::speedglm(form, data = calib, family = binomial())
		model_multivar <- glm2(form, data = calib, family = binomial())

		say('')
		print(summary(model_multivar))
		say('')

		# model did not converge
		if (!model_multivar$converged | model_multivar$boundary) {
			
			say('Model did not converge. Restarting fold...')
			timings <- timings[timings$k != k]
		
		} else { # model converged
	 
			if (demesne == 'Large') saveRDS(model_multivar, paste0(output_dir, '/model_k_', prefix(k, 2), '.rds'))
	 
			# AUC
			predict_eval <- predict(model_multivar, eval, type = 'response')
			pred_eval <- prediction(predict_eval, labels = eval$forest_loss)
			auc_obs <- performance(pred_eval, 'auc')@y.values[[1]]
			tjurs_r2_obs <- mean(predict_eval[eval$forest_loss == 1]) - mean(predict_eval[eval$forest_loss == 0])
			
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
					auc = auc_obs,
					tjurs_r2 = tjurs_r2_obs
				)
			)

			### univariate models and permutation test on multivariate model
			################################################################
			n <- nrow(eval)

			for (this_pred in predictors_selected) {

				### univariate model
				####################
				
				form <- paste0('forest_loss ~ 1 + ', this_pred)
				if (any(this_pred == 'log10_population')) form <- paste0(form, ' + I(log10_population^2)')
				if (any(this_pred == 'log10_dist_to_roads_km')) form <- paste0(form, ' + I(log10_dist_to_roads_km^2)')
				if (any(this_pred == 'log10_dist_to_rivers_km')) form <- paste0(form, ' + I(log10_dist_to_rivers_km^2)')
				if (any(this_pred == 'short_veg_density_33cells')) form <- paste0(form, ' + I(short_veg_density_33cells^2)')
				if (any(this_pred == 'forest_density_33cells')) form <- paste0(form, ' + I(forest_density_33cells^2)')

				form <- as.formula(form)
				
				# model_univar <- speedglm::speedglm(form, data = calib, family = binomial())
				model_univar <- glm2(form, data = calib, family = binomial())

				# AUC
				predict_eval <- predict(model_univar, eval, type = 'response')
				pred_eval <- prediction(predict_eval, labels = eval$forest_loss)
				auc_univar <- performance(pred_eval, 'auc')@y.values[[1]]
				tjurs_r2_univar <- mean(predict_eval[eval$forest_loss == 1]) - mean(predict_eval[eval$forest_loss == 0])
				
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
						auc = auc_univar,
						tjurs_r2 = tjurs_r2_univar
					)
				)

				### leave-out permutation test
				##############################
				
				eval_copy <- eval
				auc_permute <- tjurs_r2_permute <- rep(NA_real_, n_permute)

				# iterate permutation test
				for (i in seq_len(n_permute)) {
				
					eval_copy[ , this_pred] <- sample(eval_copy[[this_pred]], n)
					predict_eval_permute <- predict(model_multivar, eval_copy, type = 'response')

					pred_eval <- prediction(predict_eval_permute, labels = eval$forest_loss)
					auc_permute[i] <- performance(pred_eval, 'auc')@y.values[[1]]
					tjurs_r2_permute[i] <- mean(predict_eval_permute[eval$forest_loss == 1]) - mean(predict_eval_permute[eval$forest_loss == 0])
				
				}

				auc_permute <- mean(auc_permute)
				tjurs_r2_permute <- mean(tjurs_r2_permute)

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
						auc = auc_permute,
						tjurs_r2 = tjurs_r2_permute
					)
				)

			} # next predictor

			if (demesne != 'Small') {

				### univariate model: PA * country
				###################################
				
				form <- paste0('forest_loss ~ 1 + country + PA + country:PA')
				form <- as.formula(form)
				
				# model_univar <- speedglm::speedglm(form, data = calib, family = binomial())
				model_univar <- glm2(form, data = calib, family = binomial())

				# AUC
				predict_eval <- predict(model_univar, eval, type = 'response')
				pred_eval <- prediction(predict_eval, labels = eval$forest_loss)
				auc_univar <- performance(pred_eval, 'auc')@y.values[[1]]
				tjurs_r2_univar <- mean(predict_eval[eval$forest_loss == 1]) - mean(predict_eval[eval$forest_loss == 0])
				
				# get coefficient values
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
						coefficient = 'country × PA',
						coefficient_value = NA,
						auc = auc_univar,
						tjurs_r2 = tjurs_r2_univar
					)
				)

				### leave-out permutation test: country * PA
				##############################################
				
				eval_copy <- eval
				auc_permute <- tjurs_r2_permute <- rep(NA_real_, n_permute)

				# iterate permutation test
				for (i in seq_len(n_permute)) {
				
					eval_copy[ , 'country'] <- sample(eval_copy[['country']], n)
					eval_copy[ , 'PA'] <- sample(eval_copy[['PA']], n)
					predict_eval_permute <- predict(model_multivar, eval_copy, type = 'response')

					pred_eval <- prediction(predict_eval_permute, labels = eval$forest_loss)
					auc_permute[i] <- performance(pred_eval, 'auc')@y.values[[1]]
					tjurs_r2_permute[i] <- mean(predict_eval_permute[eval$forest_loss == 1]) - mean(predict_eval_permute[eval$forest_loss == 0])
				
				}

				auc_permute <- mean(auc_permute)
				tjurs_r2_permute <- mean(tjurs_r2_permute)

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
						coefficient = 'country & PA',
						coefficient_value = NA,
						auc = auc_permute,
						tjurs_r2 = tjurs_r2_permute
					)
				)


			} # not small demesne

			say('CROSS-VALIDATION FOLD ', k, ' PREDICTION RASTERS', level = 2)
			step <- paste0('Cross-validation: Make prediction raster')
				
			start <- Sys.time()
			prediction_rast <- predict(env, model_multivar, type = 'response')
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'predict()', target = 'Predictor raster', dtype = 'raster', start = start, stop = stop, k = k)
			
			start <- Sys.time()
			names(prediction_rast) <- paste0('prediction_k', prefix(k, 3))
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'names()', target = 'Predictor raster', dtype = 'raster', start = start, stop = stop, k = k)

			bigTiff <- demesne %in% c('Medium', 'Large')
			start <- Sys.time()
			writeRaster(prediction_rast, paste0(output_dir, '/prediction_k', prefix(k, 3), '.tif'), bigTiff = bigTiff, overwrite = TRUE)
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Predictor raster', dtype = 'raster', start = start, stop = stop, k = k)

			prediction_rasts[[k]] <- prediction_rast
			
			fwrite(timings, paste0(output_dir, tolower(demesne), '_fasterRaster_timings.csv'), row.names = FALSE)
			fwrite(results, paste0(output_dir, tolower(demesne), '_fasterRaster_results.csv'), row.names = FALSE)

			# # # # NB commenting this out bc mow() can't detect GSpatials that are in lists
			# # # ### CLEAN-UP DISK SPACE
			# # # #######################
			
			# # # say('REMOVE TEMPORARY FILES', level = 2)
			# # # step <- 'Remove temporary files'

			# # # start <- Sys.time()
			# # # mow('unlinked', ask = FALSE)
			# # # stop <- Sys.time()
			# # # timings <- remember(timings, step = step, fx = 'mow()', target = 'Temporary files', dtype = 'raster/vector', start = start, stop = stop)

			k <- k + 1
				
		} # model converged

	} # next fold

	### POST-PROCESS PREDICTION RASTERS
	###################################
	
	say('POST-PROCESS PREDICTION RASTERS', level = 2)
	step <- 'Post-process prediction rasters'

	# combine prediction rasters
	start <- Sys.time()
	predictions <- do.call(c, prediction_rasts)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictions))

	# average and dispersion of rasters
	start <- Sys.time()
	prediction_mean <- mean(predictions)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mean()', target = 'Prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(predictions))

	start <- Sys.time()
	prediction_sd <- stdev(predictions, pop = FALSE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'stdev()', target = 'Prediction rasters', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	prediction_cv <- prediction_sd / prediction_mean
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Mean and SD of prediction rasters', dtype = 'raster', start = start, stop = stop)

	# stack mean and CV prediction
	start <- Sys.time()
	prediction_mean_cv <- c(prediction_mean, prediction_cv)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'c()', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop)

	# name
	start <- Sys.time()
	names(prediction_mean_cv) <- c('mu', 'cv')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))

	# save prediction rasters (mean, CV)
	byLayer <- bigTiff <- demesne %in% c('Medium', 'Large')
	fn <- paste0(output_dir, tolower(demesne), '_fasterRaster_prediction_mean_cv.tif')
	start <- Sys.time()
	writeRaster(prediction_mean_cv, fn, overwrite = TRUE, byLayer = byLayer, bigTiff = bigTiff)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Mean and CV of prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))

	### CLEAN-UP DISK SPACE
	#######################
	
	say('REMOVE TEMPORARY FILES', level = 2)
	step <- 'Remove temporary files'

	keep <- list(prediction_mean, prediction_sd, prediction_cv, prediction_mean_cv)
	start <- Sys.time()
	mow('*', ask = FALSE, keep = keep)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'mow()', target = 'Temporary files', dtype = 'raster/vector', start = start, stop = stop)

	### classify prediction rasters into 9 classes based on combination of probability of forest loss and uncertainty
	#################################################################################################################
	say('CLASSIFY PREDICTION RASTERS INTO 9 CLASSES BASED ON COMBINATION OF PROBABILITY OF FOREST LOSS AND UNCERTAINTY', level = 2)
	
	step <- 'Classify predictions'
	
	# NB these steps are restricted to fasterRaster when runnig the Medium or Large demesnes because they fail on terra. As a result, a workaround must be used in terra, so the function will ot match.
	
	restricted <- if (demesne == 'Small') { NA_character_} else { 'fasterRaster' }
	
	start <- Sys.time()
	quants_mean <- global(prediction_mean, fun = 'quantile', probs = threshold_quantiles)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'global()', target = 'Mean prediction raster', dtype = 'raster', start = start, stop = stop, restricted = restricted, n = length(threshold_quantiles))

	start <- Sys.time()
	quants_cv <- global(prediction_cv, fun = 'quantile', probs = threshold_quantiles)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'global()', target = 'CV prediction rasters', dtype = 'raster', start = start, stop = stop, restricted = restricted, n = length(threshold_quantiles))

	quants_mean <- unlist(quants_mean)
	quants_cv <- unlist(quants_cv)

	# # stack mean and CV prediction
	# start <- Sys.time()
	# prediction_mean_cv <- c(prediction_mean, prediction_cv)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'c()', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop)

	# # name
	# start <- Sys.time()
	# names(prediction_mean_cv) <- c('mu', 'cv')
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'names()', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))

	# CV min/max for creating reclass table
	start <- Sys.time()
	mm_cv <- minmax(prediction_cv)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'minmax()', target = 'CV of predictions rasters', dtype = 'raster', start = start, stop = stop)

	# calculate classes based on low/medium/high probability of forest loss and uncertainty therein
	# we assign the mean predictions classes of 1 to 3, and CV clases of 10, 20, or 30, then add the two class rasters to get a final class
	mean_reclass_table <- data.frame(
		low = c(0, quants_mean[1], quants_mean[2]),
		high = c(quants_mean[1], quants_mean[2], 1),
		replace = 1:3
	)

	cv_reclass_table <- data.frame(
		low = c(0, quants_cv[1], quants_cv[2]),
		high = c(quants_cv[1], quants_cv[2], mm_cv['max', 1]),
		replace = c(10, 20, 30)
	)

	start <- Sys.time()
	prediction_mean_class <- classify(prediction_mean, mean_reclass_table, include.lowest = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'classify()', target = 'Mean prediction raster', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	prediction_cv_class <- classify(prediction_cv, cv_reclass_table, include.lowest = TRUE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'classify()', target = 'CV of predictions rasters', dtype = 'raster', start = start, stop = stop)

	start <- Sys.time()
	prediction_class <- prediction_mean_class + prediction_cv_class
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'arithmetic', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop, n = 2)

	# zone_fx <- paste0(' = if (mu <= ', quants_mean[1], ' & cv <= ', quants_cv[1], ', 0, if (mu > ', quants_mean[1], ' & mu <= ', quants_mean[2], ' & cv <= ', quants_cv[1], ', 1, if (mu > ', quants_mean[2], ' & cv <= ', quants_cv[1], ', 2, if (mu <= ', quants_mean[1], ' & cv > ', quants_cv[1], ' & cv <= ', quants_cv[2], ', 3, if (mu > ', quants_mean[1], ' & mu <= ', quants_mean[2], ' & cv > ', quants_cv[1], ' & cv <= ', quants_cv[2], ', 4, if (mu >= ', quants_mean[2], ' & cv > ', quants_cv[1], ' & cv <= ', quants_cv[2], ', 5, if (mu <= ', quants_mean[1], ' & cv > ', quants_cv[2], ', 6, if (mu > ', quants_mean[1], ' & mu <= ', quants_mean[2], ' & cv > ', quants_cv[2], ', 7, if (mu > ', quants_mean[2], ' & cv > ', quants_cv[2], ', 8, null())))))))))')

	# # NB using cores > 1 throws error: first error: 'unused arguments (quants_mean = c(0.2665024 . . .))'
	# start <- Sys.time()
	# prediction_class <- app(prediction_mean_cv, fun = zone_fx)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'app()', target = 'Mean and CV of predictions rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))
	
	# name
	start <- Sys.time()
	names(prediction_class) <- 'prediction_class'
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'names()', target = 'Predictions class raster', dtype = 'raster', start = start, stop = stop)

	# # make categorical
	# start <- Sys.time()
	# prediction_class <- as.int(prediction_class)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'as.int()', target = 'Predictions class raster', dtype = 'raster', start = start, stop = stop)

	levs <- data.frame(
		value = c(11, 12, 13, 21, 22, 23, 31, 32, 33),
		class = c(
			'low prob/low uncert',
			'medium prob/low uncert',
			'high prob/low uncert',
			'low prob/medium uncert',
			'medium prob/medium uncert',
			'high prob/medium uncert',
			'low prob/high uncert',
			'medium prob/high uncert',
			'high prob/high uncert'
		)
	)

	start <- Sys.time()
	levels(prediction_class) <- levs
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'levels()<-', target = 'Predictions class raster', dtype = 'raster', start = start, stop = stop)

	# # save prediction rasters (mean, CV)
	# byLayer <- bigTiff <- demesne %in% c('Medium', 'Large')
	# fn <- paste0(output_dir, tolower(demesne), '_fasterRaster_prediction_mean_cv.tif')
	# start <- Sys.time()
	# writeRaster(prediction_mean_cv, fn, overwrite = TRUE, byLayer = byLayer, bigTiff = bigTiff)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Mean and CV of prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))

	bigTiff <- demesne %in% c('Medium', 'Large')
	fn <- paste0(output_dir, tolower(demesne), '_fasterRaster_prediction_class.tif')
	start <- Sys.time()
	writeRaster(prediction_class, fn, overwrite = TRUE, bigTiff = bigTiff)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Class of prediction rasters', dtype = 'raster', start = start, stop = stop)

	fwrite(timings, paste0(output_dir, tolower(demesne), '_fasterRaster_timings.csv'), row.names = FALSE)

say('DONE!', deco = '^', level = 1)
say(date())
sink()
