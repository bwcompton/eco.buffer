'write.tiff' <- function(x, template, name, what = '', verbose = FALSE, timing = FALSE, type = NULL) {

   # write.tiff
   # write matrix as a geoTIFF, which bells and whistles
   # Arguments:
   #   x          matrix to write
   #   template   raster object that x conforms to (e.g., input raster)
   #   name       path and name of result
   #   what       what to call it in chatter (start with lowercase)
   #   verbose    whether to chatter
   #   timing     whether to report timing (only if chatter = TRUE)
   #   type       one of int, float, or NULL to figure it out. Alternatively, you can use one
   #              of the types specified for writeRaster
   # B. Compton, 29 May 2021
   # 12 Jan 2024: use terra, add type
   # 16 Jan 2024: use assessType, copied from rasterPrep, to get proper nodata value



   chatter(verbose, 'Writing ', what, ' grid...')
   t <- proc.time()[3]

   if(is.null(type))                                        # resolve type
      type <- ifelse(isTRUE(all.equal(x, floor(x))), 'int', 'float')
   type <- switch(type,
          int = 'INT2S',
          float = 'FLT4S',
          type
   )

   template[] <- as.vector(x)
   if(length(grep('.tif$', name, ignore.case = TRUE)) == 0) name <- paste0(name, '.tif')
   writeRaster(template, name, datatype = type, overwrite = TRUE, NAflag = assessType(type)$noDataValue)

   what <- paste(toupper(substring(what, 1, 1)), substring(what, 2), sep = '')
   chatter(verbose, what, ' grid written to ', name)

   chatter(timing & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')
}
