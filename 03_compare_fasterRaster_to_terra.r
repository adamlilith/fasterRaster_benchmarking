### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script compares the results of the same workflow performed with fasterRaster and terra.
###
### source('C:/Ecology/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/03_compare_fasterRaster_to_terra.r')
### source('E:/Adam/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/03_compare_fasterRaster_to_terra.r')
###
### CONTENTS ###
### settings ###
### compare fasterRaster and terra runtimes ###

################
### settings ###
################

	rm(list=ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	setwd(paste0(drive, '/Research/fasterRaster - Streamlined GIS in R through GRASS'))

	library(data.table) # data frames
	library(enmSdmX) # SDMing
	library(ggplot2) # graphics
	library(omnibus) # utilities
	library(terra) # GIS

say('###############################################')
say('### compare fasterRaster and terra runtimes ###')
say('###############################################')

	demesne <- 'Small'

	terra <- fread(paste0('./outputs/completed/', tolower(demesne), '_terra_timings.csv'))
	fr <- fread(paste0('./outputs/completed/', tolower(demesne), '_fasterRaster_timings.csv'))

	terra[ , ':=' (time_s = runtime_s / n)]
	fr[ , ':=' (time_s = runtime_s / n)]

	fx <- apply(fr[ , c('fx', 'datatype')], MARGIN = 1, paste, collapse = ' ')
	fr[ , method_match := ifelse(fx == 'fast() raster', 'rast()', ifelse(fx == 'vect() vector', 'vect()', terra$fx))]

	no_matches <- fr$fx[fr$method_match %notin% terra$fx]
