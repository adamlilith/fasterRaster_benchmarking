### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script conducts benchmark tests for the fasterRaster, terra, and sf packages for R.
###
### source('C:/Ecology/Drive/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/00_constants.r')
### source('E:/Adam/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/00_constants.r')
###
### CONSTANTS ###

# buffer width around Mekong basin in meters
study_region_buffer_size_m <- 40000

# names of basins for largest extent
primary_basins_large <- c('Mekong', 'Salween', 'Irrawaddy', 'Chao Phraya', 'Sittang', 'Hong (Red River)', 'Xun Jiang')

# names of countries in study region
focal_country_names <- c('China', 'India', 'Myanmar', 'Thailand', 'Laos', 'Cambodia', 'Vietnam')

# number of sites in model calibration set and in evaluation set (each), as proportion of total cells in forest in study region
# cross_valid_prop <- 0.1
cross_valid_prop <- 0.01

# number of cross-validation folds
# n_folds <- 5
n_folds <- 2

# number of times to conduct permutation test for each variable
n_permute <- 100

# thresholds (quantiles) to define low/high probability of deforestation and low/high confidence
threshold_quantiles <- c(0.25, 0.75)

### function to record timings
# df		data frame
# step		description of step
# fx	 	name of function (like: terra::vect() or fasterRaster::vect())
# target	name of target (like: "Mekong countries")
# dtype		"vector" or "raster"
# start,stop	difftime object (time elapsed)
# n 		Number of items processed (i.e., number of rasters ingested)
remember <- function(df, step, fx, target, dtype, start, stop, n = 1) {

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
			runtime_s = runtime_s
		)
	)

	omnibus::say(step, ': ', target, ': ', fx, ': ', round(runtime_s, 3), ' sec')

	rownames(df) <- NULL
	df

}
