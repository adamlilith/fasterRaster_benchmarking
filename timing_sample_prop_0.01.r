n_points_inflated
[1] 12394871
> start <- Sys.time()
> sample_sites <- spatSample(forest_2000, n_points_inflated, values = TRUE, as.points = TRUE, verbose = TRUE)
 Ingesting points...

  |                                    
  |                              |   0%
  |                                    
  |                              |   2%
  |                                    
  |=                             |   3%
  |                                    
  |=                             |   5%
  |                                    
  |==                            |   6%
  |                                    
  |==                            |   8%
  |                                    
  |===                           |  10%
  |                                    
  |===                           |  11%
  |                                    
  |====                          |  13%
  |                                    
  |====                          |  15%
  |                                    
  |=====                         |  16%
  |                                    
  |=====                         |  18%
  |                                    
  |======                        |  19%
  |                                    
  |======                        |  21%
  |                                    
  |=======                       |  23%
  |                                    
  |=======                       |  24%
  |                                    
  |========                      |  26%
  |                                    
  |========                      |  27%
  |                                    
  |=========                     |  29%
  |                                    
  |=========                     |  31%
  |                                    
  |==========                    |  32%
  |                                    
  |==========                    |  34%
  |                                    
  |===========                   |  35%
  |                                    
  |===========                   |  37%
  |                                    
  |============                  |  39%
  |                                    
  |============                  |  40%
  |                                    
  |=============                 |  42%
  |                                    
  |=============                 |  44%
  |                                    
  |==============                |  45%
  |                                    
  |==============                |  47%
  |                                    
  |===============               |  48%
  |                                    
  |===============               |  50%
  |                                    
  |===============               |  52%
  |                                    
  |================              |  53%
  |                                    
  |================              |  55%
  |                                    
  |=================             |  56%
  |                                    
  |=================             |  58%
  |                                    
  |==================            |  60%
  |                                    
  |==================            |  61%
  |                                    
  |===================           |  63%
  |                                    
  |===================           |  65%
  |                                    
  |====================          |  66%
  |                                    
  |====================          |  68%
  |                                    
  |=====================         |  69%
  |                                    
  |=====================         |  71%
  |                                    
  |======================        |  73%
  |                                    
  |======================        |  74%
  |                                    
  |=======================       |  76%
  |                                    
  |=======================       |  77%
  |                                    
  |========================      |  79%
  |                                    
  |========================      |  81%
  |                                    
  |=========================     |  82%
  |                                    
  |=========================     |  84%
  |                                    
  |==========================    |  85%
  |                                    
  |==========================    |  87%
  |                                    
  |===========================   |  89%
  |                                    
  |===========================   |  90%
  |                                    
  |============================  |  92%
  |                                    
  |============================  |  94%
  |                                    
  |============================= |  95%
  |                                    
  |============================= |  97%
  |                                    
  |==============================|  98%
  |                                    
  |==============================| 100%
 Collapsing points...

  |                                    
  |                              |   0%
  |                                    
  |====                          |  14%
  |                                    
  |=========                     |  29%
  |                                    
  |=============                 |  43%
  |                                    
  |=================             |  57%
  |                                    
  |=====================         |  71%
  |                                    
  |==========================    |  86%
  |                                    
  |==============================| 100%
 Extracting values...
 Creating GVector...
> stop <- Sys.time()
> timings <- remember(timings, step = step, fx = 'spatSample()', target = 'Forest loss raster masked to focal basin', dtype = 'raster', start = start, stop = stop)
 Define study region: Forest loss raster masked to focal basin: spatSample(): 17669.285 sec | Mon Jun 17 05:48:31 2024
> 
> # remove points in NA cells
> # NB we could use fasterRaster::complete.cases() on the GVector, but this combines two functions we want to time
> start <- Sys.time()
> values_at_points <- as.data.frame(sample_sites)
> stop <- Sys.time()
> timings <- remember(timings, step = step, fx = 'as.data.frame()', target = 'Candidate calibration/evaluation sites', dtype = 'vector', start = start, stop = stop)
 Define study region: Candidate calibration/evaluation sites: as.data.frame(): 0.023 sec | Mon Jun 17 05:48:31 2024
> 
> completes <- which(complete.cases(values_at_points))
> start <- Sys.time()
> sample_sites <- sample_sites[completes]
> stop <- Sys.time()
> timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)
 Define study region: Candidate calibration/evaluation points: subset_single_bracket: 2978.706 sec | Mon Jun 17 06:38:10 2024
> 
> # keep desired number of point
> start <- Sys.time()
> n_sample_sites <- nrow(sample_sites)
> stop <- Sys.time()
> timings <- remember(timings, step = step, fx = 'nrow()', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)
 Define study region: Candidate calibration/evaluation points: nrow(): 0.01 sec | Mon Jun 17 06:38:10 2024
> 
> keeps <- sample(n_sample_sites, cross_valid_n)
> 
> start <- Sys.time()
> sample_sites <- sample_sites[keeps]
> stop <- Sys.time()
> timings <- remember(timings, step = step, fx = 'subset_single_bracket', target = 'Candidate calibration/evaluation points', dtype = 'vector', start = start, stop = stop)
 Define study region: Candidate calibration/evaluation points: subset_single_bracket: 1877.702 sec | Mon Jun 17 07:09:29 2024
> 
