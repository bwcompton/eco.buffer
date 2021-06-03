'broccoli.kernel' <- function(threshold, kern, id, streams, flow, accumulation, cellsize) {

   # broccoli.kernel
   # Build broccoli kernels for eco.buffer: where kernels fall in stream cells
   # with flow accumulation <= threshold, expand the kernel to fill the entire
   # watershed.
   #
   # Arguments:
   #   threshold     threshold (km^2)
   #   kern          kernel value
   #   id             id id
   #   streams       stream centerlines (=1)
   #   flow          d8 flow grid
   #   accumulation  d8 flow accumulation grid
   #   cellsize      resolution of raster cells in map units (usually m)
   # Results:
   #   kern and id, extended to fill small watersheds
   # Note: downflow takes 12 min for all of MA. It takes 13 sec for 4300 selected cells - this should be fine.
   # B. Compton, 30-31 May 2021



   # print('Here we are, in broccoliland!')
   # threshold<<-threshold;kern<<-kern;id<<-id;streams<<-streams;flow<<-flow;accumulation<<-accumulation;cellsize<<-cellsize;return()

   thresh <- threshold * 1e6 / prod(cellsize)
   streams[is.na(streams)] <- 0
   f <- kern
   g <- id
   n <- x <- (id != 0) & (accumulation <= thresh) & (streams == 1)# targets: kernel, below threshold, in stream centerline
   f <- kern * x
   g <- id * x
   while(any(n)) {                                                # While we're still finding stuff,
      if(chatter)
         cat('.')
      f <- pmax(f, downflow(n * f, flow))                         #    find cells upstream from target cells and assign id ids
      g <- pmax(g, downflow(n * g, flow))                         #    do it again for kernels
      n <- (f != 0) & (x == 0)   & (accumulation <= thresh)       #    mark new cells (these are the first outside of existing kernel)
      x <- f
   }
   if(chatter)
      cat('\n')
   list(f, g)     # return kernel and seed id
}
