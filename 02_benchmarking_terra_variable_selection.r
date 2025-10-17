### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, fasterRaster, and sf packages for R.
###
### source('C:/fasterRaster_benchmarking_distributed/fasterRaster_benchmarking/02_benchmarking_terra_variable_selection.r')
###
### CONTENTS ###
### settings ###
### benchmark ###

################
### settings ###
################

	rm(list=ls())

	drive <- 'C:/'

	.libPaths(paste0(drive, '/fasterRaster_benchmarking_distributed/libraries'))
	setwd(paste0(drive, '/fasterRaster_benchmarking_distributed'))
	source('./fasterRaster_benchmarking/00_constants.r')

	memory <- 64 * 1024
	cores <- 4
	threads <- TRUE

	grass_dir <- 'C:/Program Files/GRASS GIS 8.4/'
	faster(grassDir = grass_dir, memory = memory, cores = cores, useDataTable = TRUE, verbose = TRUE)

	memmax <- 64

	terraOptions(
		memmax = memmax,
		progress = 0
	)

say('#################')
say('### benchmark ###')
say('#################')

	### options
	###########

	demesne <- 'Large'
	# demesne <- 'Medium'
	# demesne <- 'Small'

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
		
		variable_selection_inflated_n <- NA

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

		variable_selection_inflated_n <- NA

	} else if (demesne == 'Large') {

		basins_primary <- c('Mekong', 'Salween', 'Irrawaddy', 'Chao Phraya', 'Sittang')
		basins_secondary <- NA
		country_names <- c('China', 'India', 'Myanmar', 'Thailand', 'Laos', 'Cambodia', 'Vietnam')
		n_folds <- n_folds_large
		glad_lab_tile_names <- c('30N_090E', '20N_100E', '20N_090E', '10N_100E', '40N_090E', '30N_100E')
		pop_tile_names <- c('GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R5_C27', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C27', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R6_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C27', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R7_C29', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R8_C28', 'GHS_POP_E2000_GLOBE_R2023A_54009_100_V1_0_R8_C29')
		extent <- c(73, 135, 8, 54)
		aggregate_factor <- large_agg_factor
		cross_valid_n <- cross_valid_n_large

		variable_selection_inflated_n <- 998194456

	}

	# output_dir <- paste0(substr(drive, 1, 2), '/!scratch/terra_outputs_', tolower(demesne), '/')
	output_dir <- paste0('./terra_outputs_', tolower(demesne), '/')

	sink(paste0(output_dir, tolower(demesne), '_terra_benchmark_variable_selection.txt'), split = TRUE)
	say('BENCHMARK TERRA VARIABLE SELECTION ONLY')
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
	say('Assessing terra!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', post = 2)
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

	# study region
	say('Study region...')
	study_region <- vect(paste0(output_dir, '/', tolower(demesne), '_terra_study_region.gpkg'))

	# predictors
	say('Response and predictor rasters...')
	env <- rast(paste0(output_dir, '/', tolower(demesne), '_scaled_terra_response_predictors.tif'))

	### VARIABLE SELECTION
	######################

	say('VARIABLE SELECTION', level = 2)
	step <- 'Variable selection'
	
	# need to select points in sets for large demesne bc "bad_alloc" error if we select all at once
	# NB tried selecting a smaller set and rbind()'ing it to the previous set of sets, but this also yielded bad_alloc error, so now just saving to disk an rbind()'ing after the loop
	if (demesne == 'Large') {

		dirCreate(output_dir, '/point_samples')
		sets <- 1000
		for (set in 1:sets) {

			n <- if (set < sets) {
				floor(variable_selection_inflated_n / sets)
			} else {
				floor(variable_selection_inflated_n / sets) + (variable_selection_inflated_n - sets * floor(variable_selection_inflated_n / sets))
			}

			start <- Sys.time()
			this_variable_selection_points <- spatSample(study_region, n)
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)
			
			start <- Sys.time()
			writeVector(this_variable_selection_points, paste0(output_dir, '/point_samples/variable_selection_set_', prefix(set, 5), '.gpkg'), overwrite = TRUE)
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)
			
			rm(this_variable_selection_points)
			
			start <- Sys.time()
			tmpFiles(current = TRUE, orphan = FALSE, remove = TRUE)
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'tmpFiles()', target = 'Temporary files', dtype = 'raster & vector', start = start, stop = stop, restricted = 'terra')
			
			gc()
			
		}
		
	} else {
	
		start <- Sys.time()
		variable_selection_points <- spatSample(study_region, variable_selection_inflated_n)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop)
	
	}
	
	if (demesne == 'Large') {

		for (set in 1:sets) {
		
			if (set == 1) {
			
				start <- Sys.time()
				variable_selection_points <- vect(paste0(output_dir, '/point_samples/variable_selection_set_', prefix(set, 5), '.gpkg'))
				stop <- Sys.time()
				timings <- remember(timings, step = step, fx = 'vect()', target = 'Variable selection point samples', dtype = 'vector', start = start, stop = stop, restricted = 'terra')

			} else {
			
				start <- Sys.time()
				this_variable_selection_points <- vect(paste0(output_dir, '/point_samples/variable_selection_set_', prefix(set, 5), '.gpkg'))
				stop <- Sys.time()
				timings <- remember(timings, step = step, fx = 'vect()', target = 'Variable selection point samples', dtype = 'vector', start = start, stop = stop, restricted = 'terra')

				start <- Sys.time()
				variable_selection_points <- rbind(variable_selection_points, this_variable_selection_points)
				stop <- Sys.time()
				timings <- remember(timings, step = step, fx = 'rbind()', target = 'Variable selection point samples', dtype = 'vector', start = start, stop = stop, restricted = 'terra')

			}
			
			start <- Sys.time()
			tmpFiles(current = TRUE, orphan = FALSE, remove = TRUE)
			stop <- Sys.time()
			timings <- remember(timings, step = step, fx = 'tmpFiles()', target = 'Temporary files', dtype = 'raster & vector', start = start, stop = stop)
			
			gc()
			
		}
			
	}
		
	start <- Sys.time()
	variable_selection_values <- extract(env, variable_selection_points, ID = FALSE)
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'extract()', target = 'Response and predictor rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(env))

	variable_selection_values <- as.data.table(variable_selection_values)
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

	write.csv(timings, paste0(output_dir, tolower(demesne), '_terra_timings_variable_selection.csv'), row.names = FALSE)

say('DONE!', deco = '^', level = 1)
say(date())
sink()
