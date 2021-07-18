'apply.mult' <- function(x, mult, minmax) {

   # apply.mult
   # Apply multiplier to resistance matrix or vector
   # result ranges from 1 (assuming minimum in x is 0-1) to approximately mult[1] * max
   # use negative mult[1] to invert values
   # min and max come from (1) additional elements of mult, (2) minmax, or (3) from x
   # Arguments:
   #   x       resistance matrix or vector
   #   mult    [1] multipler (negative to invert values), [2] optional user-supplied min,
   #           [3] optional user-supplied max
   #   minmax  minimum and maximum of resist (from grid before clipping)
   # B. Compton, 25-26 May 2021
   # 18 Jul 2021: it was possible to get resistances < 1 with negative resistance and nonzero minimum



   mm <- c(min(x), max(x))    # min and max from data
   if(!is.null(minmax))       # better: supplied min and max
      mm <- minmax
   if(length(mult) >= 2)      # even better: from mult (user-supplied)
      mm[1] <- mult[2]
   if(length(mult) >= 3)
      mm[2] <- mult[3]

   if(mult[1] < 0)            # if inverting,
      x <- mm[2] - x

   (x - min(1, mm[1] * (mult > 0))) * abs(mult[1]) + 1
}
