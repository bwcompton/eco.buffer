'downflow' <- function(x, flow, include = flow != 0) {
   # Give value in x of cells immediately downstream from each cell in flow
   # Arguments:
   #   x        value matrix
   #   flow     D8 flow matrix
   #   include  optional matrix marking target cells. The 3x3 window around each target
   #            cell will be run to allow marking the upflow cells.
   # Result:
   #   matrix of downstream x values for each cell. Cells flowing out of the
   #     matrix get NA, as do cells without flow.
   # B. Compton, 17-18 May 2021
   
   
   
   x <- rbind(NA, cbind(NA, x, NA), NA)         # pad to prevent disaster
   z <- x * NA
   
   p <- matrix(1, 3, 3)                         # use four-neighbor rule for include
   include[is.na(include)] <- 0                 # exclude missing too
   include <- as.matrix(focal(raster(include), p, fun = 'sum', pad = TRUE, padValue = 0)) != 0      # take focal of include matrix
   
   for(i in 2:(dim(x)[1] - 1)) 
      for(j in 2:(dim(x)[2] - 1))
         if(include[i - 1, j - 1])
            z[i, j] <- switch(match(flow[i - 1, j - 1], 2 ^ (0:7), nomatch = 9), 
                              x[i, j + 1], x[i + 1, j + 1], x[i + 1, j], 
                              x[i + 1, j - 1], x[i, j - 1], x[i - 1, j - 1],
                              x[i - 1, j], x[i - 1, j + 1], NA)
   
   z[2:(dim(z)[1] - 1), 2:(dim(z)[2] - 1)]      # drop padding
}