'add.path' <- function(path, file) {

    # add.path
    # Adds path to a filename if it doesn't already have a path
    # Filename already has a path if it starts with / or <drive letter>:
    # Arguments:
    #   path    path to prepend (optionally ending with /)
    #   file    filename
    # B. Compton, 2 Apr 2021


    if(!is.null(path) & !is.null(file))
        if((nchar(file) > 0)) {
            s <- ifelse(substr(path, nchar(path), nchar(path)) != '/', '/', '')
            if(1 != length(grep('^/|^[a-zA-Z]:', file)))
                file <- paste(path, s, file, sep = '')
        }
    file
}
