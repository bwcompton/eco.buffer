'find.edges' <- function(x) {

    # find.edges
    # Find edge cells of polygons in raster x
    # B. Compton, 11 Apr 2021


    p <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)                     # use four-neighbor rule
    (focal(!is.na(x), p, pad = TRUE, padvalue = 0) <= 3) * !is.na(x)    # focal sum
}
