start <- Sys.time()
sample_sites <- spatSample(forest_2000, n_points_inflated, values = TRUE, as.points = TRUE, verbose = TRUE)
stop <- Sys.time()
timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Forest loss raster masked to focal basin', dtype = 'raster', start = start, stop = stop)

# remove points in NA cells
# NB we could use fasterRaster::complete.cases() on the GVector, but this combines two functions we want to time
start <- Sys.time()
values_at_points <- as.data.frame(sample_sites)
stop <- Sys.time()
timings <- remember(timings, step = step, fx = 'as.data.frame()', target = 'Candidate calibration/evaluation sites', dtype = 'vector', start = start, stop = stop)

completes <- which(complete.cases(values_at_points))
start <- Sys.time()
sample_sites <- sample_sites[completes]
stop <- Sys.time()
timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)

# keep desired number of point
start <- Sys.time()
n_sample_sites <- nrow(sample_sites)
stop <- Sys.time()
timings <- remember(timings, step = step, fx = 'nrow()', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)

keeps <- sample(n_sample_sites, cross_valid_n)

start <- Sys.time()
sample_sites <- sample_sites[keeps]
stop <- Sys.time()
timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)


for 0.01 prop, takes at LEAST 30 hr

