## fasterRaster BENCHMARKING
###
### This script conducts benchmark tests for the fasterRaster, fasterRaster, and sf packages for R. It differs from 02_benchmarking_fasterRaster_loop.r in that that code gathers all preceeding cross-validation rasters and conducts the post-modeling analysis, too. This script can be run independent of other folds being done elewhere.
###
### source('C:/fasterRaster_benchmarking_distributed/fasterRaster_benchmarking/02_benchmarking_fasterRaster_loop_in_chunks.r')
###
### CONTENTS ###
### settings ###
### benchmark ###

################
### settings ###
################

	rm(list=ls())

	# folds to do
	
	# k_set 7:50 # 
	
	# k_set <- 1:50
	# k_set <- 11:14 # DONE
	# k_set <- 15:19 # IP
	# k_set <- 19
	# k_set <- 20:24 # IP
	# k_set <- 24
	# k_set <- 25:30 # IP
	k_set <- 29
	# k_set <- 31:34 # DONE
	# k_set <- 35:39 # DONE
	# k_set <- 40:44 # DONE
	# k_set <- 46:50 # DONE
	
	all_ok <- function() {
		message(paste0('Running set k ', paste(k_set, collapse = ', ')))
		x <- readline('is this OK? (y/n) ')
		if (x != 'y') stop('All your base are belong to us.')
	}
	
	all_ok()
	
	drive <- 'C:/'

	.libPaths(paste0(drive, '/fasterRaster_benchmarking_distributed/libraries'))
	setwd(paste0(drive, '/fasterRaster_benchmarking_distributed'))
	source('./fasterRaster_benchmarking/00_constants.r')

	memory <- 64 * 1024
	cores <- 4

	grass_dir <- 'C:/Program Files/GRASS GIS 8.4/'
	faster(grassDir = grass_dir, memory = memory, cores = cores, useDataTable = TRUE, verbose = TRUE)

	threads <- TRUE
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

	all_ok <- function() {
		message(paste0('Demesne ', demesne))
		x <- readline('is this OK? (y/n) ')
		if (x != 'y') stop('All your base are belong to us.')
	}

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
		continuous_predictors_selected <- c('elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'short_veg_density_33cells', 'log10_population')

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
		
		continuous_predictors_selected <- c('elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'short_veg_density_33cells', 'log10_population')

		# predictors sected in initial run
		continuous_predictors_selected <- c('forest_density_33cells', 'elev_scale_11_m', 'slope_scale_7', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'log10_population')

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

		continuous_predictors_selected <- c('elev_scale_11_m', 'slope_scale_11', 'log10_dist_to_rivers_km', 'log10_dist_to_roads_km', 'short_veg_density_33cells', 'log10_population')

	}
	
	predictors_selected <- c(continuous_predictors_selected, 'forest_frag_class', 'country', 'PA')

	# output_dir <- paste0(substr(drive, 1, 2), '/!scratch/terra_outputs_', tolower(demesne), '/')
	output_dir <- paste0('./fasterRaster_outputs_', tolower(demesne), '_agg_factor_', aggregate_factor, '/')
	dirCreate(output_dir)

	k_set_nice <- paste0(k_set[1], '_to_', k_set[length(k_set)])
	sink(paste0(output_dir, tolower(demesne), '_fasterRaster_benchmark_k_set_', k_set_nice, '.txt'), split = TRUE)
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
	say('FOLD SET: ........................................................ ', paste(k_set, collapse = ', '))
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

	### LOAD PRIOR PREDICTION RASTERS, STUDY REGION, RESPONSE, and PREDICTORS
	#########################################################################

	# study region
	say('Study region...')
	study_region <- fast(paste0(output_dir, '/', tolower(demesne), '_fasterRaster_study_region.gpkg'))
	
	# predictors
	say('Response and predictor rasters...')
	forest_loss_focal_basin <- fast(paste0(output_dir, tolower(demesne), '_scaled_fasterRaster_response_predictors_forest_loss.tif'))
	names(forest_loss_focal_basin) <- 'forest_loss'
	env <- forest_loss_focal_basin

	for (pred in predictors_selected) {
	
		say(pred)
	
		this_pred <- fast(paste0(output_dir, '/', tolower(demesne), '_scaled_fasterRaster_response_predictors_', pred, '.tif'))
		names(this_pred) <- pred
		
		env <- c(env, this_pred)
	
	}
	
	### CROSS-VALIDATION SETUP
	##########################
	
	say('CROSS-VALIDATION SETUP', level = 2)
	step <- 'Cross-validation setup'

	# calculate number of sites for calibration and evaluation
	# calculate sample size... selecting more than needed bc sites can be placed in NA cells

	not_na <- not.na(forest_loss_focal_basin)

	# NB in fasterRaster, we could use nonnacell() to get this, but the underlying algorithm is the same
	non_na_cells <- global(not_na, 'sum')
	non_na_cells <- non_na_cells$sum

	cell_res <- res(forest_loss_focal_basin)
	non_na_cells_area_m2 <- prod(100000 * cell_res) * non_na_cells

	study_region_area_m2 <- expanse(study_region)

	# number of random points to draw
	# inflating by 1.2 to ensure we have enough points on non-NA cells
	# also inflating by proportion of NA cells bc we will need to discard points that fall in these cells

	# fold sample size
	inflation <- 1.2 * study_region_area_m2 / non_na_cells_area_m2
	cross_valid_inflated_n <- ceiling(inflation * cross_valid_n)

	say('non_na_cells ................... ', non_na_cells)
	say('for each fold, initially selecting ', cross_valid_inflated_n ,' points, which will be subset to ', cross_valid_n, ' after removing NAs')

	vars_selected <- c('forest_loss', continuous_predictors_selected, 'forest_frag_class', 'country', 'PA')
	predictors_selected <- c(continuous_predictors_selected, 'forest_frag_class', 'country', 'PA')

	env <- env[[vars_selected]]

	### CROSS VALIDATION
	####################

	# NB We do not select random points from the full predictor stack first. Rather, we select random points from the forest loss/no loss raster, then use it to remove points that are outside the study region (NA cells) and to select the desired number of cells of each class (forest loss/persistence). *Then* we extract from the full response/predictor stack. This is much faster than selecting points and extracting values using the full stack initially.

	# store results from modeling in this data frame
	results <- data.table()

	for (k in k_set) {

		say('CROSS-VALIDATION FOLD ', k, ' - BENCHMAKING IN CHUNKS', level = 2)
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

		# initial values
		if (demesne == 'Small') {
			n_terms <- 18
		} else if (demesne == 'Medium') {
			n_terms <- 21
		} else if (demesne == 'Large') {
			n_terms <- 26
		}
		start <- rep(0, n_terms)
		
		model_multivar <- glm2(form, data = calib, family = binomial(), start = start)

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

			# prediction_rasts[[k]] <- prediction_rast
			
			fwrite(timings, paste0(output_dir, tolower(demesne), '_fasterRaster_timings_k_', k, '.csv'), row.names = FALSE)
			fwrite(results, paste0(output_dir, tolower(demesne), '_fasterRaster_results_k_', k, '.csv'), row.names = FALSE)

			# # # # NB commenting this out bc mow() can't detect GSpatials that are in lists
			# # # ### CLEAN-UP DISK SPACE
			# # # #######################
			
			# # # say('REMOVE TEMPORARY FILES', level = 2)
			# # # step <- 'Remove temporary files'

			# # # start <- Sys.time()
			# # # mow('unlinked', ask = FALSE)
			# # # stop <- Sys.time()
			# # # timings <- remember(timings, step = step, fx = 'mow()', target = 'Temporary files', dtype = 'raster/vector', start = start, stop = stop)

		} # model converged

	} # next fold

	# # # fwrite(timings, paste0(output_dir, tolower(demesne), '_fasterRaster_timings_k_set_', k_set_nice, '.csv'), row.names = FALSE)

say('DONE!', deco = '^', level = 1)
say(date())
sink()
