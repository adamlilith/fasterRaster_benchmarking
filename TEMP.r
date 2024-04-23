library(tictoc)



cats <- sample(11291, round(0.56 * 11291)) # works up to 0.56 then fails for larger!!!
# cats <- sort(cats)
# cats <- seqToSQL(cats)
# cats <- as.character(cats)
cats <- paste(cats, collapse=',')

tic()
if (.exists('zero')) .rm('zero')
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'zero', cats = cats, flags = c('overwrite', 'verbose'))
.exists('zero')
.vectInfo('zero')
toc()

tic()
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'one', cats = '1-1000', flags = 'overwrite')
toc()

tic()
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'two', cats = '1-1000', type = 'area', flags = 'overwrite')
toc()

tic()
say('FAILS')
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'three', cats = '1-1000', type = 'boundary', flags = 'overwrite')
toc()

tic()
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'four', cats = paste(1:1000, collapse=','), flags = 'overwrite')
toc()

tic()
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'five', cats = paste(1:1000, collapse=','), type = 'area', flags = 'overwrite')
toc()

tic()
say('FAILS')
src <- .makeSourceName("v_extract", "vector")
rgrass::execGRASS("v.extract", input = sources(se_asia_basins), output = 'six', cats = paste(1:1000, collapse=','), type = 'boundary', flags = 'overwrite')
toc()

