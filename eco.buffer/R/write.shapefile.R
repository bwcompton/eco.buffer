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



   chatter(verbose, 'Writing ', what, ' shapefile...')
   t <- proc.time()[3]

   shapefile(x, file = name, overwrite = TRUE)   # and write shapefile

   what <- paste(toupper(substring(what, 1, 1)), substring(what, 2), sep = '')

   if(verbose)
      chatter(verbose, 'Result shapefile written to ', name)

   chatter(timing & verbose, '  Elapsed time = ', proc.time()[3] - t, ' s')
}
