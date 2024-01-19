'write.shapefile' <- function(x, name, what = '', verbose = FALSE, timing = FALSE) {

   # write.shapefile
   # write shapefile object, which bells and whistles
   # Arguments:
   #   x          shapefile object to write
   #   name       path and name of result
   #   what       what to call it in chatter (start with lowercase)
   #   verbose    whether to chatter
   #   timing     whether to report timing (only if chatter = TRUE)
   # B. Compton, 29 May 2021
   # 12 Jan 2024: use sf
   # 19 Jan 2024: Use backslashes in file paths



   chatter(verbose, 'Writing ', what, ' shapefile...')
   t <- proc.time()[3]

   if(length(grep('.shp$', name, ignore.case = TRUE)) == 0) name <- paste0(name, '.shp')
   suppressWarnings(st_write(x, name, delete_dsn = TRUE, quiet = TRUE))            # and write shapefile

   what <- paste(toupper(substring(what, 1, 1)), substring(what, 2), sep = '')

   if(verbose)
      chatter(verbose, 'Result shapefile written to ', gsub('/+', '\\\\', name))

   chatter(timing & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')
}
