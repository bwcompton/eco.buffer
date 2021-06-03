'downflow' <- function(value, flow) {
   # Give value of cells immediately downstream from each cell in flow
   # Arguments:
   #   value    value matrix
   #   flow     D8 flow matrix
   # Result:
   #   matrix of values downstream from each cell
   # B. Compton, 31 May 2021



   f <- rbind(0, cbind(0, flow, 0), 0)                # pad to prevent disaster
   f[is.na(f)] <- 0
   v <- as.vector(t(rbind(0, cbind(0,  value, 0), 0)))# final indices
   z <- v * 0

   c <- matrix(c(1, 1, 0, 2, 1, 1, 4, 0, 1, 8, -1, 1, 16, -1, 0, 32, -1, -1, 64, 0, -1, 128, 1, -1), 8, 3, byrow = TRUE) # control matrix: flowdir, x offset, y offset
   d <- c[, 3] * dim(f)[2] + c[, 2]                   # offset for each flow direction

   for(i in 1:8) {                                    # For each flow direction,
      j <- t(f == c[i, 1])                            #    index into this flow direction
      z[j] <- v[(1:length(v))[j] + d[i]]              #    get cell index for downstream cells in this direction
   }
   z <- matrix(z, dim(f)[1], dim(f)[2], byrow = TRUE) # reshape matrix
   z <- z[2:(dim(z)[1] - 1), 2:(dim(z)[2] - 1)]       # drop padding
   z
}
