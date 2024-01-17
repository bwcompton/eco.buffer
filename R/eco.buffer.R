#' Ecological Buffer Tool
#'
#' Delineates ecologically-defined buffers around target conservation areas using a resistant kernel
#'
#' @param seeds a shapefile of points, lines, or polygons designating conservation
#' targets
#' @param bandwidth maximum distance (m) to spread through cells of resistance = 1.
#' @param seedid name of numeric ID column in seeds shapefile.
#' @param result name of result polygon shapefile of ecological buffers.
#' @param resultgrid name of result kernel grid (optional). If resultgrid is
#' supplied, a graduated kernel representation. of buffers will be created, useful
#' for exploring parameterization and visualizing kernel-creation process.
#' @param resist landscape resistance grid. Optional if landcover and resist.table
#' are used to assign resistance by landcover type. Resistance values should range
#' from 1 to infinity. Use resist.mult to invert positive grids.
#' @param resist.table optional tab-delimited text file. Defines resistance value
#' for each landcover class (ranging from 1 to infinity). Columns must include class
#' (numeric landcover class) and resist (resistance value); additional columns may
#' include landcover name for reference (highly recommended) and comments. Resistance
#' ranges from 1 to infinity. If both resist and resist.table are supplied, the
#' values from the resist raster will be used for values that are absent in the
#' table.
#' @param default.resist value to use for classes not included in table, or NULL to
#' throw error if any values in landcover are missing from table. Nodata cells
#' (outside of the landscape) always get a high resistance.
#' @param landcover landcover grid. May be used to designate resistance values with
#' resist.table; also used with screen option.
#' @param barrier grid with resistance values for barriers, used as a complement to
#' resist or resist.table, the maximum resistance of resist/resist.table will be used
#' for each cell. Use to supply values for aquatic barriers (bridges, culverts, and
#' dams), with nodata for other cells.
#' @param passage grid with resistance values for cells that provide passage, reducing
#' resistance, used as a complement to resist or resist.table. The minimum of
#' resist/resist.table and passage will be used for each cell. Use to supply terrestrial
#' passage over or under roads (bridges and wildlife passage structures), with nodata
#' for other cells.
#' @param clip a polygon shapefile to clip the analysis to. Use this when developing
#' parameters and testing to speed runs up immensely.
#' @param fullextent if TRUE, produces a result grid at the full extent of the
#' landscape. This runs more slowly, so use the default of FALSE unless you have a
#' good reason not to. If TRUE, all grids are clipped to the reference grid
#' (landcover if it exists, otherwise, reference grid) to enforce alignment.
#' @param resist.mult multiplier on resistance grid, such that such that resistance
#' values that originally range from 0 to n or 1 to n will range from 1 to approximately
#' n * resist.mult after multiplying. Use a negative multiplier to invert the resistances
#' when using a positive grid, where higher values denote lower resistance. The minimum
#' and maximum resistances are needed for this procedure. They are obtained from the entire
#' raster before any clipping. If you are running on tiled data, or data that otherwise
#' don't represent the entire range of resistance values, you need to supply the overall
#' minimum and maximum as the 2nd and 3rd elements. The default multiplier is 1, thus
#' resistance values are used directly (or 1 is added if the minimum is 0).
#' @param resist.table.mult multiplier on resistances from a table. See resist.mult.
#' @param barrier.mult multiplier on barrier resistances. See resist.mult.
#' @param passage.mult multiplier on passage resistances. See resist.mult.
#' @param save.resist specify a result TIFF to write the realized resistance grid to,
#' for assessing complex combinations of resistance sources.
#' @param path base path prepended to input and result names/paths. Inputs that
#' include a complete path (starting with / or a drive letter) don't use path. This
#' option helps keep the inputs cleaner, and makes for easy switching to different
#' sets of inputs and results.
#' @param simplify if TRUE, simplify result polygons; if FALSE, polygons exactly
#' match raster.
#' @param simplify.tolerance polygon smoothing parameter used if simplify = TRUE.
#' Larger values give simpler polygons.
#' @param density build kernels for every nth cell to speed things up. By default,
#' density = 1, and resistant kernels are built from each edge cell in seeds; a
#' larger value will decrease runtime (by density2) at the cost of precision. You can
#' get away with higher values for density when using larger bandwidths.
#' @param expand distance to expand seeds (m). Use with screen to apply lines to wide
#' streams, for example. If expanded seeds overlap, the overlapping area will be
#' arbitrarily assigned the seedid of one of the seeds.
#' @param screen limit seeds to these landcover classes (requires landcover). Use
#' this if seeds are sloppy, e.g. when designating stream cores from vector data that
#' don't correspond exactly to streams in landcover, or to exclude development
#' classes from conservation target polygons.
#' @param broccoli include entire watershed above point in stream if watershed area
#' at point is <= x km2 (requires flow, accumulation, and streams grids). Used to
#' include the entirety of small watersheds.
#' @param streams stream centerline grid, used with the broccoli option.
#' @param flow flow direction grid, used with broccoli option. Grid is a standard D8
#' grid.
#' @param accumulation flow accumulation grid, used with broccoli option.
#' @param verbose set to FALSE to suppress informational chatter.
#' @param timing display overall timing message if TRUE, and also intermediate timing
#' messages if 'all' (verbose must be TRUE too).
#'
#' @details
#' The ecological buffer tool **eco.buffer** is a stand-alone R package for delineating
#' ecologically-defined buffers around target conservation areas. Targeted areas are
#' designated by point, line, or polygon shapefiles. Buffers are based on resistant
#' kernels with flexible parameters to accommodate terrestrial or aquatic settings.
#' Landscape resistance can be defined by a landcover raster and a table of classes
#' and resistance values, directly from resistance rasters, or from a combination of
#' a resistance table and rasters. Results are a polygon shapefile of ecological
#' buffers and an optional geoTIFF raster of the resistant kernels.
#'
#' The buffer tool is based on resistant kernels (Compton et. al 2007), which have
#' been used in a number of conservation applications since 2003,
#' including estimating local and regional connectivity (McGarigal et al. 2018) and
#' building terrestrial and aquatic conservation cores in Designing Sustainable
#' Landscapes/Nature's Network (McGarigal et al. 2017); they have also been used in
#' TNC's Resilient Sites for Terrestrial Conservation and Massachusetts Natural
#' Heritage's Living Waters and BioMap 2.
#'
#' @section Notes:
#' 1. Resistance values must range from 1 to infinity. The spread value starts in
#' each focal cell (all edge cells of each seed) at bandwidth / cell size. At each
#' cell, the cell's resistance x multiplier is subtracted from the spread value. For
#' example, a bandwidth of 5000 m when the cell size is 30 m gives a spread value of
#' 166.67. The spread will stop once it has passed through cells with a cumulative
#' resistance * multiplier of 166.67. Resistances greater than or equal to this value
#' will stop the spread at a single cell, thus these cells act as complete barriers.
#' 2. Raster inputs may be either Arc grids or geoTIFFs (other formats will likely
#' work).
#' 3. The seeds shapefile may be singlepart or multipart.
#' 4. The seed id field (specified with the seedid option) will be preserved in the
#' resulting buffer polygon. When result buffers overlap, the id of the seed with a
#' shortest cost-distance to each point will be used.
#' 5. When using the CAPS landcover for terrestrial cores, make sure to set the
#' resistance of classes 60 (Bridge or culvert) and 61 (Dam) to the maximum road
#' resistance, as these classes interrupt roads, thus with a lower resistance, they
#' can allow a kernel to spread across roads where they occur.
#' 6. When using the CAPS landcover for aquatic cores, you may want to use the
#' Aquatic Barriers (abarriers) grid as a secondary resistance grid using the barriers
#' option, to assign resistance to each bridge or culvert and each dam based on estimates
#' of their aquatic passability.
#' 7. If both a resistance raster (resist option) and resistance table (resist.table and
#' landcover options) are supplied, the table is used for all classes in the table, and
#' the raster is used for any classes not in the table (default.resist will be ignored).
#' 8. When testing and developing parameters, runtime will be much faster if you
#' limit seeds to a relatively small geographic area, or use the clip option select a
#' small area of the landscape.
#' 9. All input grids will be snapped and clipped to the reference grid--the landcover,
#' if supplied, or else the resistance grid.
#' 10. When using the broccoli option with a flow accumulation grid in terms of cells
#' (as is the CAPS grid), convert to km^2 with 1e6 / cellsize ^ 2. For a 30 m grid, use
#' flow accumulation / 1111 to get km^2.
#' 11. Note that when using polygons as seeds, tiny or skinny polygons that cover
#' less than half of a cell may not be properly captured when converting seeds to
#' raster. eco.buffer ensures that at least the centroid (moved inside of the
#' polygon) will be captured, but long polygons that are narrower than a cell may be
#' represented by only a single cell or a disjunct series of cells. It's always good
#' practice to use the resultgrid option to save a raster version of the result, and
#' look at cells with a value of 1.0 to see where the seed polygons ended up in the
#' raster representation.
#'
#' @section References:
#' Compton, B.W., K. McGarigal, S.A. Cushman, and L.R. Gamble. 2007. A resistant-kernel
#' model of connectivity for amphibians that breed in vernal pools. Conservation
#' Biology 21:788-799. \doi{10.1111/j.1523-1739.2007.00674.x}.
#'
#' McGarigal, K., B.W. Compton, E.B. Plunkett, W.V. DeLuca, J. Grand, E. Ene, and
#' S.D. Jackson. 2018. A landscape index of ecological integrity to inform landscape
#' conservation. Landscape Ecology 33:1029-1048. \doi{10.1007/s10980-018-0653-9}.
#'
#' McGarigal K., B.W. Compton, E.B. Plunkett, W.V. DeLuca, and J. Grand. 2017.
#' Designing sustainable landscapes: landscape conservation design. Report to the
#' North Atlantic Conservation Cooperative, US Fish and Wildlife Service, Northeast
#' Region.
#' \url{http://landeco.umass.edu/web/lcc/dsl/technical/DSL_documentation_landscape_design.pdf}
#' @section Author:
#' Bradley W. Compton <bcompton@@umass.edu>
#' @export
#  C++ code that does resistant kernels. Package source is at https://github.com/ethanplunkett/gridprocess
#' @import gridprocess
#  GIS processing for raster data
#' @import terra
#' @import sf
#' @importFrom utils read.table
#' @examples
#' ### Set up temporary directory for examples
#' require(eco.buffer)
#' dp <- paste(shortPathName(system.file('exampledata', package='eco.buffer')), '/.', sep = '')
#' dir <- tempdir()
#' if(!file.exists(dir)) dir.create(dir)
#' file.copy(dp, dir, recursive=TRUE)
#' cat('Example data and results will be in', dir)
#'
#' ### 1. terrestrial kernels (creates test1.shp and testg1.tif)
#' eco.buffer('seed_points', bandwidth = 2000, landcover = 'capsland.tif',
#' resist.table = 'resistance.txt', result = 'test1', resultgrid = 'testg1.tif', path = dir)
#'
#' ### 2. WMA poly example (creates WMAtest2.shp)
#' eco.buffer('wma_seeds', 5000, landcover = 'capsland.tif', resist = 'iei.tif', resist.mult = -30,
#'      resist.table = 'resist_dev.txt', result = 'WMAtest2', path = dir)
#'
#' ### 3. stream cores (creates stream_test3.tif)
#' eco.buffer('stream_seeds', bandwidth = 3000, landcover = 'capsland.tif',
#'      resist.table = 'resist_streams.txt', default.resist = 999, barrier = 'abarriers.tif',
#'      barrier.mult = 100, result = 'stream_test3', simplify = FALSE, path = dir)
# B. Compton, 2 Apr 2021-31 May 2021
# 11 Jan 2024: revise to use sf and terra



'eco.buffer' <- function(seeds, bandwidth, seedid = 'Id', result = NULL, resultgrid = NULL, resist = NULL,
                         resist.table = NULL, default.resist = NULL, landcover = '',
                         barrier = NULL, passage = NULL, resist.mult = 1, resist.table.mult = 1,
                         barrier.mult = 1, passage.mult = 1, save.resist = NULL, path = '', density = 1,
                         expand = NULL, screen = NULL, broccoli = NULL, flow = NULL, accumulation = NULL,
                         streams = NULL, clip = NULL, fullextent = FALSE, simplify = TRUE,
                         simplify.tolerance = 30, verbose = TRUE, timing = FALSE) {



   t0 <- proc.time()[3]

   if(any(timing == 'all'))
      timing <- c(TRUE, TRUE)
   else {
      if(length(timing) == 1)
         timing <- c(timing, FALSE)
   }

   # Process file names and check for required files
   if(!is.null(path))
      cat('path = ', path, '\n', sep = '')
   seeds <- check.file(seeds, path, 'Seeds shapefile', verbose, ext = 'shp', require = TRUE)
   clip <- check.file(clip, path, 'Clip shapefile', verbose, ext = 'shp')
   resist <- check.file(resist, path, 'Resistance grid', verbose)
   resist.table <- check.file(resist.table, path, 'Resistance table', verbose)
   landcover <- check.file(landcover, path, 'Landcover grid', verbose)
   barrier <- check.file(barrier, path, 'Barrier grid', verbose)
   passage <- check.file(passage, path, 'Passage grid', verbose)
   streams <- check.file(streams, path, 'Streams grid', verbose)
   flow <- check.file(flow, path, 'Flow grid', verbose)
   accumulation <- check.file(accumulation, path, 'Flow accumulation grid', verbose)
   result <- check.file(result, path, 'Result shapefile', verbose, ext = 'shp', result = TRUE)
   resultgrid <- check.file(resultgrid, path, 'Result grid', verbose, ext = 'tif', result = TRUE)
   save.resist <- check.file(save.resist, path, 'Saved resistance grid', verbose, ext = 'tif', result = TRUE)


   # Check file dependencies
   if(!is.null(resist.table) & is.null(landcover))
      stop('When supplying resist.table, you must also supply landcover')
   if(is.null(resist) & is.null(resist.table))
      stop('You must supply either resist or resist.table (or both)')
   if(is.null(result) & is.null(resultgrid) & is.null(save.resist))
      stop('You must supply result, resultgrid, or both...or at least save.result')
   if(!is.null(screen) & is.null(landcover))
      stop('When screen option is supplied, you must also supply landcover')
   if(!is.null(broccoli) & (is.null(flow) | is.null(accumulation) | is.null(streams)))
      stop('When broccoli is being used, you must supply all three of flow, accumulation, and streams')



   # ------------------------ READ AND CLIP DATA ------------------------
   chatter(verbose, '\nReading seeds shapefile...')
   seedshape <- st_read(add.path(path, seeds), quiet = TRUE)

   if(is.na(match(seedid, names(seedshape))))
      stop('Seed ID column (seedid = \'', seedid, '\') not found in shapefile ', seeds)

   if(all(st_is(seedshape, 'POINT')) & density != 1)
      stop('When seeds are points, denisty must be 1 (otherwise some seeds may be lost)')

   if(!any(is.na(suppressWarnings(as.numeric(seedshape[[seedid]])))))
      seedshape[[seedid]] <- as.numeric(seedshape[[seedid]])      # force ID to numeric
   else
      stop('The field referred to by seedid (currently ',seedid, ') must be numeric')

   # if expand option is included, buffer seeds
   if(!is.null(expand) && expand > 0) {
      q <- buffer(seedshape, expand, dissolve = FALSE)    # buffer seeds

      # I think these 3 lines are unnecessary. The just move seedid to 1st column. They're preserved in buffering,
      # and subsequent references use the name, not position, so why would I care?
      # names(q)[1] <- seedid
      # q[[seedid]] <- seedshape[[seedid]]                  # and transfer seedid
      # seedshape <- q
   }

   if(!is.null(clip)) {
      chatter(verbose, 'Reading clip shapefile...')
      clipshape <- st_read(add.path(path, clip), quiet = TRUE)
   }

   # Read resistance grid, if supplied to get reference
   if(!is.null(resist)) {
      chatter(verbose, 'Reading resistance grid...')
      rg <- rast(add.path(path, resist))
      ref <- rg
   }

   # Read landcover grid (optional) and get reference
   if(!is.null(landcover)) {
      chatter(verbose, 'Reading landcover grid...')
      land <- rast(add.path(path, landcover))
      ref <- land
   }

   # Clip grids to seeds + buffer for speed-up, unless fullextent, in which case we clip to reference grid
   t <- proc.time()[3]
   if(!is.null(clip)){               # if clip is set, use it
      chatter(verbose, 'Clipping to clip box...')
      shapes <- suppressWarnings(st_intersection(seedshape, clipshape))   # clip shapefile

      ref <- crop(ref, clipshape)
   }
   else if(!fullextent) {             # otherwise, clip to seeds + bandwidth, unless fullextent, in which case clip to full reference grid
      chatter(verbose, 'Clipping grids to match ', ifelse(fullextent, 'reference grid...', 'seeds shapefile...'))
      q <- st_bbox(seedshape)
      q <- q + c(-bandwidth, bandwidth, -bandwidth, bandwidth)
      ref <- crop(ref, q)
   }

   # clip all grids for alignment, to clip, seeds + bandwidth, or fullextent, depending
   if(!is.null(resist)) {
      resist.minmax <- minmax(rg)
      r <- as.matrix(crop(rg, ref), wide = TRUE)
   }

   if(!is.null(landcover))
      land <- as.matrix(crop(land, ref), wide = TRUE)

   # read and clip barrier grid
   if(!is.null(barrier)) {
      barrg <- rast(barrier)
      barrier.minmax <- minmax(barrg)
      barrg <- as.matrix(crop(barrg, ref), wide = TRUE)
   }

   # read and clip passage grid
   if(!is.null(passage)) {
      passg <- rast(passage)
      passage.minmax <- minmax(passg)
      passg <- as.matrix(crop(passage, ref), wide = TRUE)
   }

   # Read and clip streams, flow, and accumulation
   if(!is.null(broccoli)) {
      streamg <- as.matrix(crop(rast(streams), ref), wide = TRUE)
      flowg <- as.matrix(crop(rast(flow), ref), wide = TRUE)
      accumg <- as.matrix(crop(rast(accumulation), ref), wide = TRUE)
   }
   chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')



   # ------------------------ BUILD RESISTANCE GRID ------------------------
   chatter(verbose, 'Getting resistance values...')
   t <- proc.time()[3]

   # Get resistance from grid
   if(!is.null(resist))
      r <- rx <- apply.mult(r, resist.mult, resist.minmax)     # apply multiplier and save bonus copy to combine with table

   # Get resistances from table
   if(!is.null(resist.table)) {
      res.table <- read.table(resist.table, sep = '\t', header = TRUE)
      resist.table.minmax <- minmax(res.table$resist)
      res.table$resist <- apply.mult(res.table$resist, resist.table.mult, resist.table.minmax)     # apply multiplier
      r <- matrix(res.table$resist[match(land, res.table$class)], dim(land)[1], dim(land)[2], byrow = FALSE)
      v <- sort(unique(as.vector(land))) # look for missing classes
      v <- v[!v %in% res.table$class]

      if(length(v) != 0 )                       # if any landcover values missing from table,
         if(!is.null(resist)) {                 #    if we have a resistance grid,
            r[land %in% v] <- rx[land %in% v]   #       use grid values where we have nothing from table
         } else                                 #    else, use default resistance or give an error
            if(is.null(default.resist)) {
               stop('Classes in landcover not specified in resistance table with no default.resist: ', paste(v, collapse = ', '))
            } else
               r[land %in% v] <- default.resist   # classes not listed get default.resist
   }

   r[is.na(r)] <- 1e6         # NA cells (e.g., outside of landscape) always get a huge resistance


   # Get barrier and passage resistances
   if(!is.null(barrier))
      r <- pmax(r, apply.mult(barrg, barrier.mult, barrier.minmax), na.rm = TRUE)

   if(!is.null(passage))
      r <- pmin(r, apply.mult(passg, passage.mult, passage.minmax), na.rm = TRUE)

   chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')


   # convert vector to grid
   chatter(verbose, 'Converting shapefile to grid...')
   t <- proc.time()[3]

   z <- rasterize(seedshape, ref, field = seedid) + 1       # add 1 to preserve zero ids

   # make sure we've captured tiny polygons by converting (inside) centroids too
###   q <- SpatialPointsDataFrame(gPointOnSurface(seedshape, byid = TRUE), seedshape@data)   ### make sure we're using seedid!
### this is so much simpler. I think we're good, right? Needs testing
   q <- st_point_on_surface(seedshape)
   q <- rasterize(q, ref, field = seedid) + 1               # add 1 again
   z[] <- pmax(as.matrix(z, wide = TRUE), as.matrix(q, wide = TRUE), na.rm = TRUE)

   if(all(is.na(as.matrix(z, wide = TRUE))))
      stop('Empty seeds. Perhaps you were too heavy-handed with clip?')


   # if screen, set seed cells in screened-out landcover classes to naught
   if(!is.null(screen)) {
      z <- z * t(land) %in% screen
      z[z == 0] <- NA
   }


   # find edge cells--we'll build kernels only from these
   x <- find.edges(z)

   # convert to matrix with 0s instead of NAs. Keep z for spatial reference.
   x <- as.matrix(x, wide = TRUE)       # edge cells
   x[is.na(x)] <- 0
   w <- as.matrix(z, wide = TRUE)
   w[is.na(w)] <- 0

   chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')



   # ------------------------ BUILD RESISTANT KERNELS ------------------------
   chatter(verbose, 'Building resistant kernels...')
   t <- proc.time()[3]
   h <- ceiling(bandwidth / res(z)[1])    # bandwidth in terms of cells
   y <- as.matrix(!is.na(z)) * h          # polygon insides
   c <- seq(-h, h)                        # window index vector
   for(i in 1:dim(x)[1])                  # For each row in landscape,
      if(i %% density == 0)               #   if not skipping this row
         for(j in 1:dim(x)[2])            #    For each column,
            if((j %% density == 0) & (x[i, j] != 0)) {  #      if not skipping this column and we have a seed here,
               k <- i + c                #       take square of data to work with, possibly clipped at edges
               k <- k[e1 <- (k >= 1) & (k <= dim(x)[1])]
               l <- j + c
               l <- l[e2 <- (l >= 1) & (l <= dim(x)[2])]
               q <- rawspread(r[k, l], h, match(i, k), match(j, l))
               w[k, l] <- ifelse(q > y[k, l], (q != 0) * z[i, j, 1], w[k, l]) #       seed id
               y[k, l] <- pmax(y[k, l], q)                                    #       kernel values

            }
   y <- y / max(y)      # rescale kernels 0 to 1



   # ------------------------ BROCCOLI ------------------------
   if(!is.null(broccoli)) {  # if broccoli option, build broccoli watersheds
      chatter(verbose, 'Expanding broccoli kernels...')
      t <- proc.time()[3]
      q <- broccoli.kernel(broccoli, y, w, streamg, flowg, accumg, res(z))
      y <- q[[1]]     # kernel
      w <- q[[2]]     # seed id
      chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')
   }


   y[y == 0] <- NA    # kernel values
   w[w == 0] <- NA    # seed id
   w <- w - 1         # recover zero ids




   # ------------------------ SAVE RESULTS ------------------------
   chatter(verbose, 'Saving results...')
   if(!is.null(result)) {
      z[] <- w

      chatter(verbose, 'Raster to poly...')
      t <- proc.time()[3]
      p <- rasterToPolygons(z, dissolve = TRUE)     # convert to polygon
      chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')

      if(simplify) {
         chatter(verbose, 'Simplifying...')
         t <- proc.time()[3]
         p <- gSimplify(p, tol = simplify.tolerance)       # smooth
         chatter(timing[2] & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')
      }
      write.shapefile(p, result, 'result', verbose, timing[2])
   }

   if(!is.null(resultgrid))
      write.tiff(y, z, resultgrid, 'result kernel', verbose, timing[2])

   if(!is.null(save.resist))
      write.tiff(r, z, save.resist, 'saved resistance', verbose, timing[2])

   chatter(timing[1] & verbose, '  Total elapsed time = ', proc.time()[3] - t0, ' s')
}
