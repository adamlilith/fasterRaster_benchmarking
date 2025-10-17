# source('C:/fasterRaster_benchmarking_distributed/fasterRaster_benchmarking/TEMP.r')

	### classify prediction rasters into 9 classes based on combination of probability of forest loss and uncertainty
	#################################################################################################################
	say('CLASSIFY PREDICTION RASTERS INTO 9 CLASSES BASED ON COMBINATION OF PROBABILITY OF FOREST LOSS AND UNCERTAINTY', level = 2)

	step <- 'Classify predictions'
	
	### NB For medium-sized region, terra fails on global(..., fun = quantile) with the error: Error in h(simpleError(msg, call)) : 
	# error in evaluating the argument 'x' in selecting a method for function 'quantile': data is too long
	# So, we will subsample the raster to get quantiles.
	
	if (demesne != 'Small') {

		start <- Sys.time()
		quantile_points <- spatSample(study_region, variable_selection_inflated_n)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Study region vector', dtype = 'vector', start = start, stop = stop, restricted = 'terra')
	
		# mean
		start <- Sys.time()
		quants_mean <- extract(prediction_mean, quantile_points, ID = FALSE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'extract()', target = 'Mean prediction raster', dtype = 'raster', start = start, stop = stop, n = nlyr(env), restricted = 'terra')
		
		# cv
		start <- Sys.time()
		quants_cv <- extract(prediction_cv, quantile_points, ID = FALSE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'extract()', target = 'CV prediction raster', dtype = 'raster', start = start, stop = stop, n = nlyr(env), restricted = 'terra')
		
		quants_mean <- quantile(quants_mean[ , 'mean'], threshold_quantiles, na.rm = TRUE)
		quants_cv <- quantile(quants_cv[ , 'cv'], threshold_quantiles, na.rm = TRUE)
	
	
	} else {
		
		start <- Sys.time()
		quants_mean <- global(prediction_mean, fun = quantile, probs = threshold_quantiles, na.rm = TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'global()', target = 'Mean prediction raster', dtype = 'raster', start = start, stop = stop, n = length(threshold_quantiles))

		start <- Sys.time()
		quants_cv <- global(prediction_cv, fun = quantile, prob = threshold_quantiles, na.rm=TRUE)
		stop <- Sys.time()
		timings <- remember(timings, step = step, fx = 'global()', target = 'CV prediction rasters', dtype = 'raster', start = start, stop = stop, n = length(threshold_quantiles))
		
		quants_mean <- unlist(quants_mean)
		quants_cv <- unlist(quants_cv)

	}

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

	# # save prediction rasters (mean, sd)
	# fn <- paste0(output_dir, tolower(demesne), '_terra_prediction_mean_cv.tif')
	# start <- Sys.time()
	# prediction_mean_cv<- writeRaster(prediction_mean_cv, fn, overwrite = TRUE)
	# stop <- Sys.time()
	# timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Mean and CV of prediction rasters', dtype = 'raster', start = start, stop = stop, n = nlyr(prediction_mean_cv))

	fn <- paste0(output_dir, tolower(demesne), '_terra_prediction_class.tif')
	start <- Sys.time()
	prediction_class <- writeRaster(prediction_class, fn, overwrite = TRUE, datatype = 'INT2U')
	stop <- Sys.time()
	timings <- remember(timings, step = step, fx = 'writeRaster()', target = 'Class of prediction rasters', dtype = 'raster', start = start, stop = stop)

	fwrite(timings, paste0(output_dir, tolower(demesne), '_terra_timings_k_start_', k_start, '.csv'), row.names = FALSE)

say('DONE!', deco = '^', level = 1)
say(date())
sink()
