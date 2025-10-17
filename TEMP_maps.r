# source('C:/Kaji/Research/fasterRaster - Streamlined GIS in R through GRASS/fasterRaster_benchmarking/TEMP_maps.r')

		map <- ggplot() +
			layer_spatial(study_region, fill = 'gray80')

			map <- map +
				layer_spatial(pred, aes(fill = after_stat(band1))) +
				scale_fill_manual(
					name = 'Forest',
					values = c('Persistence' = 'forestgreen', 'Loss' = 'red'),
					na.value = NA,
					na.translate = FALSE
				) +
				ggtitle('Forest') +
				annotation_scale(location = 'bottomleft') +
				annotation_north_arrow(location = 'bottomright', which_north = 'true')

			if (demesne == 'medium') {
			
				map <- map + theme(
					legend.position = 'inside',
					legend.position.inside = c(0.5, 0.6)
				)
			
			}

			map <- map + 
				theme(
					axis.ticks = element_blank(),
					panel.background = element_blank(),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					axis.line = element_blank(),
					axis.text = element_blank()
				)

print(map)
