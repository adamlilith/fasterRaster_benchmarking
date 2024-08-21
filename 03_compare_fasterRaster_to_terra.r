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

	library(cowplot) # graphics
	library(data.table) # data frames
	library(enmSdmX) # SDMing
	library(ggplot2) # graphics
	library(omnibus) # utilities
	library(terra) # GIS

say('###############################################')
say('### compare fasterRaster and terra runtimes ###')
say('###############################################')

	demesne <- 'Small'
	# demesne <- 'Large'

	fr <- fread(paste0('./outputs/completed/fasterRaster_outputs_', tolower(demesne), '/', tolower(demesne), '_fasterRaster_timings.csv'))
	terra <- fread(paste0('./outputs/completed/terra_outputs_', tolower(demesne), '/', tolower(demesne), '_terra_timings.csv'))
	
	fr$datatype[fr$datatype == 'vector/raster'] <- 'raster/vector'
	terra$datatype[terra$datatype == 'vector/raster'] <- 'raster/vector'
	
	# terra[ , ':=' (time_s = runtime_s / n)]
	# fr[ , ':=' (time_s = runtime_s / n)]

	fx <- apply(fr[ , c('fx', 'datatype')], MARGIN = 1, paste, collapse = ' ')
	fr[ , method_match := ifelse(fx == 'fast() raster', 'rast()', ifelse(fx == 'vect() vector', 'vect()', terra$fx))]

	terra <- terra[fx != 'fillHoles()']

	# no_matches <- fr$fx[fr$method_match %notin% terra$fx]
	# stopifnot(length(no_matches) == 0)
	
	collated <- data.table(
		step_fr = fr$step,
		step_terra = terra$step,
		fx_fr = fr$fx,
		fx_terra = terra$fx,
		target_fr = fr$target,
		target_terra = terra$target,
		datatype_fr = fr$datatype,
		datatype_terra = terra$datatype,
		n_fr = fr$n,
		n_terra = terra$n,
		runtime_s_fr = fr$runtime_s,
		runtime_s_terra = terra$runtime_s
	)
	
	# collated$fx_fr[collated$fx_fr == 'subset_single_brackets'] <- '[]'
	# collated$fx_terra[collated$fx_terra == 'subset_single_brackets'] <- '[]'
	
	# collated$fx_fr[collated$fx_fr == 'subset_double_brackets'] <- '[[]]'
	# collated$fx_terra[collated$fx_terra == 'subset_double_brackets'] <- '[[]]'
	
	collated[ , ratio := runtime_s_terra / runtime_s_fr]
	collated[ , log_ratio := log10(ratio)]
	collated <- collated[order(ratio)]
	
	collated[ , fx := ifelse(fx_fr == fx_terra, fx_fr, paste0(fx_fr, '/', fx_terra))]

	# plot

	time_breaks <- c(1/3600, 1/60, 1, 60, 60 * 10, 3600, 3600 * 3, 3600 * 24)
	break_labels <- c('1/3600 s', '1/60 s', '1 s', '1 min', '10 min', '1 hr', '3 hr', '1 day')
	
	total_runtime_hr_fr <- sum(collated$runtime_s_fr) / 3600
	total_runtime_hr_terra <- sum(collated$runtime_s_terra) / 3600
	
	total_runtime_hr_fr <- round(total_runtime_hr_fr, 2)
	total_runtime_hr_terra <- round(total_runtime_hr_terra, 2)
	
	runtimes_s <- c(collated$runtime_s_terra, collated$runtime_s_fr)
	lims <- c(min(runtimes_s), max(runtimes_s))
	
	compare_all <- ggplot(collated, aes(x = runtime_s_terra, y = runtime_s_fr, col = datatype_fr, shape = datatype_fr)) +
		geom_abline(slope = 1, intercept = 0, color = 'gray') +
		geom_point() +
		scale_x_continuous(
			trans = 'log10',
			breaks = time_breaks,
			labels = break_labels
		) +
		scale_y_continuous(
			trans = 'log10',
			breaks = time_breaks,
			labels = break_labels
		) +
		coord_fixed(xlim = lims, ylim = lims) + 
		xlab('terra runtime') + ylab('fasterRaster runtime') +
		geom_text(
			label = collated$fx,
			nudge_x = 0, nudge_y = 0.25,
			hjust = 'inward',
			check_overlap = TRUE,
			show.legend = FALSE
		) +
		annotate(
			'text',
			x = 0.7 * lims[1],
			y = lims[2] - 0.35 * diff(lims),
			label = paste('Total runtime\nfasterRaster: ', total_runtime_hr_fr, ' hr\nterra: ', total_runtime_hr_terra, ' hr'),
			hjust = 0
		) + 
		ggtitle(paste(demesne, 'demesne (all functions)'))
		
	collated_long <- collated[runtime_s_fr > 60 | runtime_s_terra > 60]
	runtimes_s <- c(collated_long$runtime_s_terra, collated_long$runtime_s_fr)
	lims <- c(min(runtimes_s), max(runtimes_s))

	compare_long <- ggplot(collated_long, aes(x = runtime_s_terra, y = runtime_s_fr, col = datatype_fr, shape = datatype_fr)) +
		geom_abline(slope = 1, intercept = 0, color = 'gray') +
		geom_point() +
		scale_x_continuous(
			trans = 'log10',
			breaks = time_breaks,
			labels = break_labels
		) +
		scale_y_continuous(
			trans = 'log10',
			breaks = time_breaks,
			labels = break_labels
		) +
		coord_fixed(xlim = lims, ylim = lims) + 
		xlab('terra runtime') + ylab('fasterRaster runtime') +
		geom_text(
			label = collated_long$fx,
			nudge_x = 0, nudge_y = 0.11,
			hjust = 'inward',
			check_overlap = TRUE,
			show.legend = FALSE
		) +
		ggtitle(paste(demesne, 'demesne (functions taking >60 s for at least one package)'))
		
	compares <- plot_grid(compare_all, compare_long, nrow = 1)
	print(compares)

