'minmax' <- function(x) {


   # minmax
   # Give minimum and maximum values of raster x
   # B. Compton, 26 May 2021


   x <- as.vector(x)
   c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}
