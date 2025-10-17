### fasterRaster BENCHMARKING
### Adam B. Smith | Missouri Botanical Garden, Saint Louis Missouri, USA | adam.smith@mobot.org | 2023-01
###
### This script compares the results of the same workflow performed with fasterRaster and terra.
###
### source('C:/Kaji/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/03_compare_fasterRaster_to_terra.r')
###
### CONTENTS ###
### settings ###
### compare fasterRaster and terra runtimes ###
### map of predictions ###
### map of effect of aggregation on distance calculations ###

################
### settings ###
################

	rm(list=ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/fasterRaster - Streamlined GIS in R through GRASS'))

	library(cowplot) # graphics
	library(data.table) # data frames
	library(enmSdmX) # SDMing
	library(fasterRaster) # GIS
	library(geodata) # GIS data
	library(ggplot2) # graphics
	library(ggrepel) # labels in ggplots
	library(ggspatial) # spatial ggplot
	library(omnibus) # utilities
	library(patchwork) # graphics
	library(readxl) # open Excel documents
	library(terra) # GIS
	
# say('###############################################')
# say('### compare fasterRaster and terra runtimes ###')
# say('###############################################')

# 	# should functions that were parts of folds be averaged?
# 	# collapse_folds <- FALSE # no
# 	collapse_folds <- TRUE # yes

# 	# factor by which cells were aggregated before applying distance()
# 	# agg_factor <- 32
# 	# agg_factor <- 128
# 	# agg_factor <- 256
# 	agg_factor <- 512

# 	seed <- 3 # will affect label placement in dot plots bc geom_text() is random

# 	# demesne <- 'Small'
# 	demesne <- 'Medium'
# 	# demesne <- 'Large'

# 	### load and pre-process data
# 	#############################
	
# 	collated <- read_xlsx(paste0('./', tolower(demesne), '_combined_timings_agg_factor_', agg_factor, '.xlsx'), sheet = 'timings')
# 	collated <- as.data.table(collated)
# 	unneeded <- which(!grepl(names(collated), pattern = 'match_'))
# 	collated <- collated[ , ..unneeded]

# 	# re-assign functions to different targets
# 	collated$datatype_fr[collated$fx_fr == 'spatSample()'] <- 'vector'
# 	collated$datatype_terra[collated$fx_terra == 'spatSample()'] <- 'vector'

# 	collated$datatype_fr[collated$fx_fr == 'extract()'] <- 'raster/vector'
# 	collated$datatype_terra[collated$fx_terra == 'extract()'] <- 'raster/vector'

# 	collated$datatype_fr[collated$datatype_fr == 'vector/raster'] <- 'raster/vector'
# 	collated$datatype_terra[collated$datatype_terra == 'vector/raster'] <- 'raster/vector'
# 	collated$datatype_terra[collated$datatype_terra == 'raster & vector'] <- 'raster/vector'

# 	collated$fx_terra[collated$fx_terra == 'nrow()/ngeom()'] <- 'ngeom()/nrow()'

# 	# average runtime across same step/function of each fold
# 	if (collapse_folds) {

# 		collated_folds <- collated[!is.na(k_fr)]
		
# 		# Calculate mean runtimes across same functions used on each fold
# 		collated_folds <- collated_folds[
# 			, .(
# 				step_fr = first(step_fr),
# 				fx_fr = first(fx_fr),
# 				target_fr = first(target_fr),
# 				datatype_fr = first(datatype_fr),
# 				n_fr = mean(n_fr, na.rm = TRUE),
# 				runtime_s_fr = mean(runtime_s_fr, na.rm = TRUE),
# 				restricted_fr = first(restricted_fr),
# 				k_fr = first(k_fr),
# 				step_terra = first(step_terra),
# 				fx_terra = first(fx_terra),
# 				target_terra = first(target_terra),
# 				datatype_terra = first(datatype_terra),
# 				n_terra = mean(n_terra, na.rm = TRUE),
# 				runtime_s_terra = mean(runtime_s_terra, na.rm = TRUE),
# 				restricted_terra = first(restricted_terra),
# 				k_terra = first(k_terra)
# 			),
# 			by = .(step_fr, step_terra, fx_fr, fx_terra)
# 		]

# 		collated_names <- names(collated)
# 		collated_folds <- collated_folds[ , ..collated_names]
	
# 		collated_sans_folds <- collated[step_fr %notin% c('Cross-validation: Select and extract from calibration and evaluation points', 'Cross-validation: Make prediction raster')]

# 		post <- which(collated_sans_folds$step_fr == 'Post-process prediction rasters')[1]

# 		collated_pre <- collated_sans_folds[1:(post - 1)]
# 		collated_post <- collated_sans_folds[post:nrow(collated_sans_folds)]

# 		collated <- rbind(collated_pre, collated_folds, collated_post)

# 	}

# 	cols <- c('step_fr', 'fx_fr', 'target_fr', 'datatype_fr', 'n_fr', 'runtime_s_fr', 'restricted_fr', 'k_fr')
# 	fr <- collated[ , ..cols]

# 	cols <- c('step_terra', 'fx_terra', 'target_terra', 'datatype_terra', 'n_terra', 'runtime_s_terra', 'restricted_terra', 'k_terra')
# 	terra <- collated[ , ..cols]

# 	names(fr) <- gsub(names(fr), pattern = '_fr', replacement = '')
# 	names(terra) <- gsub(names(terra), pattern = '_terra', replacement = '')

# 	# collated$fx_fr[collated$fx_fr == 'subset_single_brackets'] <- '[]'
# 	# collated$fx_terra[collated$fx_terra == 'subset_single_brackets'] <- '[]'
	
# 	# collated$fx_fr[collated$fx_fr == 'subset_double_brackets'] <- '[[]]'
# 	# collated$fx_terra[collated$fx_terra == 'subset_double_brackets'] <- '[[]]'
	
# 	### dot plot of terra timings vs fasterRaster timings with one dot per function call
# 	####################################################################################

# 	collated[ , ratio := runtime_s_terra / runtime_s_fr]
# 	collated[ , log_ratio := log10(ratio)]
# 	collated <- collated[order(ratio)]
	
# 	collated[ , fx := ifelse(fx_fr == fx_terra, fx_fr, paste0(fx_fr, '/', fx_terra))]

# 	time_breaks <- c(1/3600, 1/60, 1, 60, 60 * 10, 3600, 3600 * 3, 3600 * 24)
# 	break_labels <- c('1/3600 s', '1/60 s', '1 s', '1 min', '10 min', '1 hr', '3 hr', '1 day')
	
# 	total_runtime_hr_fr <- sum(collated$runtime_s_fr) / 3600
# 	total_runtime_hr_terra <- sum(collated$runtime_s_terra) / 3600
	
# 	total_runtime_hr_fr <- round(total_runtime_hr_fr, 2)
# 	total_runtime_hr_terra <- round(total_runtime_hr_terra, 2)
	
# 	runtimes_s <- c(collated$runtime_s_terra, collated$runtime_s_fr)
# 	lims <- c(min(runtimes_s), max(runtimes_s))
	
# 	compare_all <- ggplot(
# 			collated,
# 			aes(x = runtime_s_terra, y = runtime_s_fr, col = datatype_fr, shape = datatype_fr)
# 		) +
# 		geom_abline(slope = 1, intercept = 0, color = 'gray') +
# 		geom_point() +
# 		guides(
# 			color = guide_legend(title = 'Object class(es)', nrow = 1),
# 			shape = guide_legend(title = 'Object class(es)', nrow = 1)
# 		) +
# 		scale_x_continuous(
# 			trans = 'log10',
# 			breaks = time_breaks,
# 			labels = break_labels
# 		) +
# 		scale_y_continuous(
# 			trans = 'log10',
# 			breaks = time_breaks,
# 			labels = break_labels
# 		) +
# 		coord_fixed(xlim = lims, ylim = lims) + 
# 		xlab('terra runtime') + ylab('fasterRaster runtime') +
# 		geom_text(
# 			label = collated$fx,
# 			nudge_x = 0, nudge_y = 0.25,
# 			hjust = 'inward',
# 			check_overlap = TRUE,
# 			show.legend = FALSE,
# 			size = 3
# 		) +
# 		# annotate(
# 		# 	'text',
# 		# 	x = 0.7 * lims[1],
# 		# 	y = lims[2] - 0.35 * diff(lims),
# 		# 	label = paste('Total runtime\nfasterRaster: ', total_runtime_hr_fr, ' hr\nterra: ', total_runtime_hr_terra, ' hr'),
# 		# 	hjust = 0
# 		# ) + 
# 		ggtitle('b) All functions') +
# 		theme(
# 			plot.title = element_text(size = 13),
# 			axis.title = element_text(size = 11),
# 			axis.text = element_text(size = 9),
# 			axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
# 		)
		
# 	### dot plot of terra timings vs fasterRaster timings with one dot per function call
# 	### for all functions with runtime > threshold for at least one package
# 	####################################################################################

# 	collated_long <- collated[runtime_s_fr > 60 | runtime_s_terra > 60]
# 	runtimes_s <- c(collated_long$runtime_s_terra, collated_long$runtime_s_fr)
# 	lims <- c(min(runtimes_s), max(runtimes_s))

# 	# condense repeated functions that cluster into single point for labeling
# 	# each point is still show in the graph
# 	if (demesne == 'Small') {

# 		# subset_single_bracket
# 		collated_long_temp <- collated_long
# 		ssb <- collated_long_temp[fx == 'subset_single_bracket']
# 		ssb <- data.table(
# 			runtime_s_fr = mean(ssb$runtime_s_fr),
# 			runtime_s_terra = mean(ssb$runtime_s_terra),
# 			datatype_fr = 'vector',
# 			fx = 'subset_single_bracket'
# 		)

# 		# spatSample
# 		collated_long_temp <- collated_long
# 		ss <- collated_long_temp[datatype_fr == 'raster' & runtime_s_fr < 3600 & runtime_s_terra < 6]
# 		ss <- data.table(
# 			runtime_s_fr = mean(ss$runtime_s_fr),
# 			runtime_s_terra = mean(ss$runtime_s_terra),
# 			datatype_fr = 'raster',
# 			fx = 'spatSample()'
# 		)

# 		max_overlaps <- 37

# 	} else if (demesne == 'Medium') {
	
# 		max_overlaps <- 95
	
# 	}

# 	compare_long <- ggplot(collated_long, aes(x = runtime_s_terra, y = runtime_s_fr, col = datatype_fr, shape = datatype_fr)) +
# 		geom_abline(slope = 1, intercept = 0, color = 'gray') +
# 		geom_point() +
# 		guides(
# 			color = guide_legend(title = 'Object class(es)'),
# 			shape = guide_legend(title = 'Object class(es)')
# 		) +
# 		scale_x_continuous(
# 			trans = 'log10',
# 			breaks = time_breaks,
# 			labels = break_labels
# 		) +
# 		scale_y_continuous(
# 			trans = 'log10',
# 			breaks = time_breaks,
# 			labels = break_labels
# 		) +
# 		coord_fixed(xlim = lims, ylim = lims) + 
# 		xlab('terra runtime') + ylab('fasterRaster runtime') +
# 		geom_text_repel(
# 			label = collated_long$fx,
# 			nudge_x = 0, nudge_y = 0.11,
# 			hjust = 'inward',
# 			force_pull = 0.1,
# 			show.legend = FALSE,
# 			size = 3,
# 			max.overlaps = max_overlaps,
# 			max.time = 1,
# 			seed = seed
# 		) +
# 		# geom_text_repel(
# 		# 	data = ssb,
# 		# 	label = ssb$fx,
# 		# 	nudge_x = 0, nudge_y = 0.11,
# 		# 	hjust = 'inward',
# 		# 	force = 5,
# 		# 	show.legend = FALSE,
# 		# 	size = 3,
# 		# 	max.overlaps = 37,
# 		# 	max.time = 1
# 		# ) +
# 		# geom_text_repel(
# 		# 	data = ss,
# 		# 	label = ss$fx,
# 		# 	nudge_x = 0, nudge_y = 0.11,
# 		# 	hjust = 'inward',
# 		# 	force = 5,
# 		# 	show.legend = FALSE,
# 		# 	size = 3,
# 		# 	max.overlaps = 37,
# 		# 	max.time = 1
# 		# ) +
# 		ggtitle('c) Functions taking >1 min for at least one package') +
# 		theme(
# 			plot.title = element_text(size = 12),
# 			axis.title = element_text(size = 11),
# 			axis.text = element_text(size = 9),
# 			axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
# 		)
		
# 	compares <- plot_grid(
# 		compare_all + theme(legend.position = 'none'),
# 		compare_long + theme(legend.position = 'none'),
# 		nrow = 2,
# 		align = 'v'
# 	)

# 	### stacked bar plot of total runtime, one bar per package
# 	##########################################################

# 	stacked_terra <- terra[ , c('fx', 'datatype', 'runtime_s')]
# 	stacked_fr <- fr[ , c('fx', 'datatype', 'runtime_s')]

# 	stacked_terra <- stacked_terra[order(stacked_terra$runtime_s)]
# 	stacked_fr <- stacked_fr[order(stacked_fr$runtime_s)]

# 	stacked_terra <- stacked_terra[order(datatype)]
# 	stacked_fr <- stacked_fr[order(datatype)]

# 	stacked_terra[ , runtime_hr := runtime_s / 3600]
# 	stacked_fr[ , runtime_hr := runtime_s / 3600]

# 	stacked_terra[ , cumulative_runtime_hr := cumsum(runtime_hr)]
# 	stacked_fr[ , cumulative_runtime_hr := cumsum(runtime_hr)]

# 	stacked_terra[ , midpoint_hr := NA_real_]
# 	stacked_fr[ , midpoint_hr := NA_real_]
# 	for (i in 2:nrow(stacked_terra)) {
	
# 		stacked_terra$midpoint_hr[i] <-
# 			stacked_terra$cumulative_runtime_hr[i - 1] + 0.5 * stacked_terra$runtime_hr[i]
	
# 		stacked_fr$midpoint_hr[i] <-
# 			stacked_fr$cumulative_runtime_hr[i - 1] + 0.5 * stacked_fr$runtime_hr[i]
	
# 	}

# 	stacked_terra[ , package := 'terra']
# 	stacked_fr[ , package := 'fasterRaster']

# 	stacked <- rbind(stacked_terra, stacked_fr)

# 	stacked$fx_thinned <- stacked$fx
# 	if (demesne == 'Small') {
# 		if (!collapse_folds) {
# 			stacked$fx_thinned[stacked$runtime_hr < 1/4] <- ''
# 		} else {
# 			stacked$fx_thinned[stacked$runtime_hr < 0.05] <- ''
# 		}
# 		label_size <- 3
# 	} else if (demesne == 'Medium') {
# 		if (!collapse_folds) {
# 			stacked$fx_thinned[stacked$runtime_hr < 4.5] <- ''
# 		} else {
# 			stacked$fx_thinned[stacked$runtime_hr < 0.5] <- ''
# 		}
# 		label_size <- 2.8
# 	}

# 	stacked$datatype <- factor(stacked$datatype, levels = rev(c('raster', 'raster/vector', 'vector')))

# 	totals <- aggregate(runtime_hr ~ package, stacked, sum)
# 	totals$runtime_hr_nice <- paste(round(totals$runtime_hr, 2), 'hr')

# 	if (demesne == 'Small') {
# 		xlim <- c(0, max(totals$runtime_hr))
# 		# xlim[2] <- if (!collapse_folds) {
# 		# 	roundTo(xlim[2], 2, ceiling)
# 		# } else {
# 		# 	roundTo(xlim[2], 1, ceiling)
# 		# }
# 	} else if (demesne == 'Medium') {
# 		xlim <- c(0, max(totals$runtime_hr))
# 		# xlim[2] <- roundTo(1.02 * xlim[2], 2, ceiling)	
# 	}

# 	# total <- ggplot(stacked, aes(x = package, y = runtime_hr, fill = datatype)) +
# 	total <- ggplot(stacked, aes(x = package, y = runtime_hr, fill = datatype, color = datatype)) +
# 		geom_col(
# 			position = 'stack',
# 			# color = 'gray40',
# 			linewidth = 0.25
# 		) +
# 		scale_fill_manual(
# 			# values = c('raster' = '#66c2a5', 'vector' = '#fc8d62', 'raster/vector' = '#8da0cb')
# 			values = c('raster' = alpha('#F8766D', 0.7), 'vector' = alpha('#619CFF', 0.7), 'raster/vector' = alpha('#00BA38', 0.7))
# 		) +
# 		scale_color_manual(
# 			# values = c('raster' = '#1b9e77', 'vector' = '#d95f02', 'raster/vector' = '#7570b3')
# 			# values = c('raster' = 'darkred', 'vector' = 'darkblue', 'raster/vector' = 'forestgreen')
# 			values = c('raster' = '#F8766D', 'vector' = '#619CFF', 'raster/vector' = '#00BA38')
# 		) +
# 		geom_text(
# 			data = stacked,
# 			mapping = aes(
# 				x = package,
# 				y = midpoint_hr
# 			),
# 			label = stacked$fx_thinned,
# 			hjust = 0.5, vjust = 0.5,
# 			show.legend = FALSE,
# 			color = 'black',
# 			size = label_size,
# 			angle = 0
# 			# angle = 90
# 		) +
# 		annotate(
# 			geom = 'text',
# 			x = totals$package, y = totals$runtime_hr,
# 			label = totals$runtime_hr_nice, 
# 			# hjust = -0.2,
# 			vjust = -0.5,
# 			size = 3.5
# 		) +
# 		labs(x = NULL, y = 'Runtime (hr)') +
# 		# coord_flip() +
# 		ylim(xlim[1], xlim[2]) +
# 		ggtitle(paste0('a) Total runtime for\n     ', tolower(demesne), ' focal region\n     with runtimes\n     ', ifelse(collapse_folds, '', 'not '), 'averaged\n     across folds')) +
# 		theme(
# 			legend.position = 'none',
# 			plot.title = element_text(size = 13),
# 			axis.title = element_text(size = 12),
# 			axis.ticks.y = element_blank(),
# 			panel.background = element_blank(),
# 			panel.grid.major = element_blank(),
# 			panel.grid.minor = element_blank(),
# 			axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0.5, size = 10),
# 			axis.line.x = element_line(color = 'black'),
# 			axis.text.x = element_text(size = 10)
# 		)

# 	### mutual legend
# 	#################

# 	legend <- get_legend(
# 		# create some space to the left of the legend
# 		compare_all + theme(legend.box.margin = margin(0, 0, 0, 0))
# 	)
		
# 	compares_legend <- plot_grid(compares, legend, rel_heights = c(12, 1), ncol = 1)
# 	# together <- plot_grid(total, compares, ncol = 1, rel_heights = 0.1, 0.4)

# 	# together <- total / (compares_legend) +
# 		# plot_layout(heights = c(1.4, 4))

# 	together <- plot_grid(total, compares_legend, rel_widths = c(4, 8), ncol = 2)

# 	# together <- total + (compares_legend) +
# 		# plot_layout(widths = c(4, 8))

# 	ggsave(together, filename = paste0('./outputs_fasterRaster/runtime_comparison_', tolower(demesne), '_agg_factor_', agg_factor, ifelse(collapse_folds, '_folds_collapsed', ''), '.png'), width = 7, height = 10.5, dpi = 600, bg = 'white')

say('##########################')
say('### map of predictions ###')
say('##########################')

	### are we just writing code or wanting to make the "official" images at highest resolution?
	# prelim <- TRUE # writing code
	prelim <- FALSE # high resolution

	if (prelim) say('Using preliminary settings--plots will be low resolution!!!!!')

	### color of areas in study region that are NA in the rasters
	study_region_col <- 'gray80'

	study_region_small <- vect('./outputs_fasterRaster/fasterRaster_outputs_small_agg_factor_128/small_fasterRaster_study_region.gpkg')
	study_region_medium <- vect('./outputs_fasterRaster/fasterRaster_outputs_medium_agg_factor_512/medium_fasterRaster_study_region.gpkg')
	study_region_large <- vect('./outputs_fasterRaster/fasterRaster_outputs_large_agg_factor_1/large_fasterRaster_study_region.gpkg')
	
	# pred_small <- rast('./outputs_fasterRaster/fasterRaster_outputs_small_agg_factor_128/small_fasterRaster_prediction_mean_cv_mu.tif')
	# pred_medium <- rast('./outputs_fasterRaster/fasterRaster_outputs_medium_agg_factor_512/medium_fasterRaster_prediction_mean_cv_mu.tif')
	pred_large <- rast('./outputs_fasterRaster/fasterRaster_outputs_large_agg_factor_1/large_fasterRaster_prediction_mean_cv_mu.tif')

	# pred_small <- pred_small[[1]]

	countries <- c('India', 'Bangladesh', 'China', 'Cambodia', 'Laos', 'Thailand', 'Vietnam', 'Myanmar', 'Indonesia', 'Nepal', 'Bhutan')
	se_asia <- gadm(countries, level = 0, resolution = 2, path = 'C:/!scratch/')
	study_region_large <- study_region_large - se_asia[se_asia$COUNTRY == 'India']

	width <- 3 * 800
	height <- 3 * 900
	maxcell <- if (prelim) { 1E6 } else { 8E6 }

	png(paste0('./outputs_fasterRaster/map_probability_of_loss_2000_large.png'), res = 600, width = width, height = height)
	
		cols <- colorRampPalette(c('#fff5f0', '#fcbba1', '#fb6a4a', '#cb181d', '#67000d'))(100)

		par(oma = rep(0, 4), mar = rep(0, 4), bg = NA, col = NA)
		plot(study_region_large, col = study_region_col, border = 'blue', legend = FALSE, lwd = 0.2, background = NA)
		plot(se_asia, lwd = 0.2, col = 'gray95', add = TRUE)
		plot(study_region_large, border = 'black', col = study_region_col, add = TRUE)
		plot(pred_large, maxcell = maxcell, col = cols, legend = FALSE, range = c(0, 1), add = TRUE)
		plot(se_asia, lwd = 0.2, col = NA, add = TRUE)
		plot(study_region_large, border = 'black', add = TRUE)

	dev.off()

	png(paste0('./outputs_fasterRaster/map_nested_study_regions.png'), res = 600, width = width, height = height)

		par(oma = rep(0, 4), mar = rep(0, 4), bg = NA, col = NA)
		plot(study_region_large, border = 'blue', axes = FALSE, background = NA)
		usr <- par('usr')
		extent <- ext(usr)
		extent <- as.polygons(extent, crs = crs(se_asia))
		plot(extent, col = 'white', lwd = 1, add = TRUE)
		plot(se_asia, lwd = 0.2, col = 'gray90', add = TRUE)
		plot(study_region_large, col = '#9ecae1', border = NA, add = TRUE)
		plot(study_region_medium, col = '#4292c6', border = NA, add = TRUE)
		plot(study_region_small, col = '#08306b', border = NA, add = TRUE)
		plot(se_asia, lwd = 0.8, col = NA, add = TRUE)
		plot(extent, lwd = 1, add = TRUE)

	dev.off()

	# legend
	png(paste0('./outputs_fasterRaster/map_legend.png'), res = 600, width = width, height = height)
		par(bg = NA, cex = 0.4)
		plot(pred_large, col = cols, mar = c(0, 0, 0, 4), axes = FALSE, background = NA, box = NA)
	dev.off()

# say('########################')
# say('### tutorial example ###')
# say('########################')

# 	# Tutorial example for the article. Main output is maps/graphs of the inputs and outputs.

# 	# gr_dir <- '/Applications/GRASS-8.4.app/Contents/Resources'
# 	gr_dir <- 'C:/Program Files/GRASS GIS 8.4'
# 	# gr_dir <- '/usr/local/grass'

# 	faster(grassDir = gr_dir)

# 	madElev <- fastData("madElev")
# 	madCover <- fastData("madCover")
# 	madRivers <- fastData("madRivers")

# 	elev <- fast(madElev)
# 	cover <- fast(madCover)
# 	rivers <- fast(madRivers)

# 	river_buff <- buffer(rivers, 1000)
# 	elev_mask <- mask(elev, river_buff)
# 	geomorphs <- geomorphons(elev_mask)
# 	# plot(geomorphs, main = "Geomorphons")
# 	geomorph_freqs <- freq(geomorphs)
# 	geomorph_freqs

# 	# We use ggplot2 and ggspatial for making maps for the article.

# 	elevSR <- rast(elev)
# 	coverSR <- rast(cover)
# 	riversSV <- vect(rivers)

# 	dirCreate('./tutorials')

# 	# elevation and rivers
# 	png('./tutorials/elev_rivers.png', res = 600, width = 1200, height = 1400)
	
# 	par(mfrow = c(1, 2))
	
# 	plot(elev, main = '')
# 	plot(rivers, col = 'blue', add = TRUE)

# 	plot(geomorphs, main = '')

# 	dev.off()

# 	# geomorphon frequency
# 	geomorph_freqs <- geomorph_freqs[order(count), ]
# 	geomorph_freqs$geomorphon <- factor(geomorph_freqs$geomorphon, levels = geomorph_freqs$geomorphon)

# 	gms <- ggplot(geomorph_freqs, aes(x = geomorphon, y = count, fill = geomorphon)) +
# 		geom_bar(stat = 'identity') +
# 		labs(x = 'Geomorphon', y = 'Frequency') +
# 		theme_minimal(base_size = 14) +
# 		theme(
# 			legend.position = 'none',
# 			axis.text.x = element_text(angle = 45, hjust = 1)
# 		)

# 	ggsave(gms, filename = './tutorials/geomorphon_frequencies_barchart.png', width = 10, height = 5, dpi = 600, bg = 'white')

# say('#############################################################')
# say('### map of effect of aggregation on distance calculations ###')
# say('#############################################################')

# 	nonagged <- rast('./outputs_fasterRaster/fasterRaster_outputs_small_agg_factor_1/small_fasterRaster_response_predictors_log10_dist_to_roads_km.tif')
# 	agged <- rast('./outputs_fasterRaster/fasterRaster_outputs_small_agg_factor_128/small_fasterRaster_response_predictors_log10_dist_to_roads_km.tif')
# 	diff <- nonagged - agged

# 	nonagged <- ggplot() +
# 		layer_spatial(nonagged) +
# 		ggtitle('Not aggregated') +
# 		guides(fill = guide_colorbar(title = 'log10(distance to roads, km)'))

# 	agged <- ggplot() +
# 		layer_spatial(agged) +
# 		ggtitle('Aggregated by factor of 128') +
# 		guides(fill = guide_colorbar(title = 'log10(distance to roads, km)'))

# 	diff <- ggplot() +
# 		layer_spatial(diff) +
# 		ggtitle('Not aggregated - aggregated') +
# 		guides(fill = guide_colorbar(title = 'Difference'))

# 	maps <- nonagged + agged + diff
# 	ggsave(maps, filename = './outputs_fasterRaster/effect_of_aggregation_on_distance_to_roads.tif', dpi = 300, width = 18, height = 8)

say('FINIS', deco = '~', level = 1)
