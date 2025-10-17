# fasterRaster_benchmarking
This repository contains code for benchmarking <a href = 'https://github.com/adamlilith/fasterRaster'>`fasterRaster`</a>. The `fasterRaster` and `terra` workflows are nearly identical and are intended to be run independently of one another.  

`00_constants.r`: Constants and functions used in subsequent scripts.  
`01_preparing_data.r`: Prepares data for benchmarking.  
`02_benchmarking_fasterRaster.r`: Entire workflow based on `fasterRaster`.  
`02_benchmarking_fasterRaster_loop.r`: Only the model loop and post-loop analysis part of the workflow.  
`02_benchmarking_fasterRaster_loop_in_chunks.r`: Code to run crossvalidation subsets of the model loop.  
`02_benchmarking_terra.r`:   Entire workflow based on `terra`.  
`02_benchmarking_terra_loop.r`: Only the model loop and post-loop analysis part of the workflow.  
`03_compare_fasterRaster_to_terra.r`: Post-workflow analyses comparing the two workflows.  

~
