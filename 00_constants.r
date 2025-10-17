### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, terra, and sf packages for R.
###
### source('C:/Ecology/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/00_constants.r')
### source('C:/Adam/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/00_constants.r')
###
### SETUP ###
### CONSTANTS ###
### FUNCTIONS ###

#############
### SETUP ###
#############

	library(data.table) # fast data frames
	library(enmSdmX) # SDMing
	library(omnibus) # utilities
	library(ROCR) # AUC
	library(glm2) # GLM
	library(terra) # GIS
	library(usdm) # predictor selection
	library(glm2) # GLM with better fitting
	library(omnibus) # utilities

	# devtools::load_all(paste0(drive, '/Kaji/R/fasterRaster'))
	library(fasterRaster) # rasters and vectors

#################
### CONSTANTS ###
#################

	# buffer width around Mekong basin in meters
	study_region_buffer_size_m <- 60000

	# # names of basins for largest extent
	# primary_basins_large <- c('Mekong', 'Salween', 'Irrawaddy', 'Chao Phraya', 'Sittang', 'Hong (Red River)', 'Xun Jiang')

	# # names of countries in study region
	# focal_country_names <- c('China', 'India', 'Myanmar', 'Thailand', 'Laos', 'Cambodia', 'Vietnam')

	# number of sites in each of model calibration set and evaluation set
	# cross_valid_prop_small_fasterRaster <- 0.05
	# cross_valid_prop_medium_fasterRaster <- 0.01
	# cross_valid_prop_large_fasterRaster <- 0.001

	cross_valid_n_small <- 200000
	cross_valid_n_medium <- 800000
	cross_valid_n_large <- 3200000

	# number of cross-validation folds
	# NB 2025-02-25 Sometimes speedglm does not converge. I am declariong that I will thus use the 40 first models that do converge, discarding the rest.
	n_folds_small <- 50
	n_folds_medium <- 50
	n_folds_large <- 50

	# factor by which to aggregate cells for calculation of distance metrics
	# small_agg_factor <- 32 # too long for terra on small demesne
	# small_agg_factor <- 128
	small_agg_factor <- 128 * 2^2
	# medium_agg_factor <- 512 # started at aggregation factor of 64, but this was too slow for terra
	medium_agg_factor <- 512 * 2^2
	# large_agg_factor <- 2048 # 256 --> >2 weeks for terra
	large_agg_factor <- 1 # for fasterRaster for "real" analysis only

	# take this times number of cross-validation points to use for variable selection
	inflation_for_variable_selection <- 2

	# number of times to conduct permutation test for each variable
	n_permute <- 100

	# thresholds (quantiles) to define low/high probability of deforestation and low/high confidence
	threshold_quantiles <- c(0.25, 0.75)

#################
### FUNCTIONS ###
#################

	### function to record timings
	# df		data frame
	# step		description of step
	# fx	 	name of function (like: terra::vect() or fasterRaster::vect())
	# target	name of target (like: "Mekong countries")
	# dtype		"vector" or "raster"
	# start,stop	difftime object (time elapsed)
	# n 		Number of items processed (i.e., number of rasters ingested)
	# restricted Either NA or 'terra' or 'fasterRaster'... indicates if this step is only used for testing the respective package
	# k			Either NA_integer_ (default) or the k-fold number
	remember <- function(df, step, fx, target, dtype, start, stop, n = 1, restricted = NA_character_, k = NA_integer_) {

		if (is.null(df)) df <- data.frame()

		runtime <- stop - start
		units <- attr(runtime, 'units')
		runtime <- as.vector(runtime)
		
		scalar <- if (units == 'secs') {
			1
		} else if (units == 'mins') {
			60
		} else if (units == 'hours') {
			3600
		} else if (units == 'days') {
			24 * 3600
		} else if (units == 'weeks') {
			7 * 24 * 3600
		}
		
		runtime_s <- scalar * runtime
		
		df <- rbind(
			df,
			data.frame(
				step = step,
				fx = fx,
				target = target,
				datatype = dtype,
				n = n,
				runtime_s = runtime_s,
				restricted = restricted,
				k = k
			)
		)

		omnibus::say(step, ': ', target, ': ', fx, ': ', round(runtime_s, 3), ' sec | ', date())

		rownames(df) <- NULL
		df

	}

