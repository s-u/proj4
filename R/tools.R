# Copyright (c) 2007 Simon Urbanek
#
# proj4 R package, License: GPL v2

# convert projection argument into PROJ.4 arguments string
# accepts either a named vector c(proj='merc',units='m'),
# a vector of parameters c('+proj=merc','+units=m')
# or a single string "+proj= +units=m"
.proj2char <- function(proj) {
    if (length(names(proj))) {
        proj <- paste('+',names(proj),'=',proj,sep='',collapse='\n')
        # remove spaces in all arguments
        proj <- gsub('\n',' ',gsub(' ','',proj))
    } else {
        if (length(proj) > 1)
            proj <- paste(as.character(proj), collapse=' ')
    }
    if (!is.character(proj)) proj <- as.character(proj)
    proj
}
