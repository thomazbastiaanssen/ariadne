#' Import modules from a file or a database
#' 
#' \code{importModules} retrieves modules information from a file or database.
#' 
#' @param map.file \code{Character vector}. One or more paths to custom mapping
#'   files or one or more of the available databases (\code{"ChocoPhlAn"}).
#'
#' 
#' @name importModules
NULL

ModuleDatabases <- list(
    GMM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GMMs.v1.07.txt",
    GBM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GBMs.v1.0.txt"
)

#' @rdname importModules
#' @export
setMethod("importModules", signature = c(module.file = "character"),
    function(module.file, merge = TRUE, message = TRUE){
        if( !is.logical(merge) ){
            stop("'merge' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !is.logical(message) ){
            stop("'message' must be TRUE or FALSE.", call. = FALSE)
        }
        modules <- lapply(module.file, .import_modules, message = message)
        # Merge mapping files
        if( merge ){
            modules <- unlist(modules, recursive = FALSE)
        }
        return(modules)
    }
)

.import_modules <- function(module.file, message){
    # Whether to use package or custom modules
    if( module.file %in% names(ModuleDatabases) ){
        # Cache database
        module.file <- .getCache(ModuleDatabases[[module.file]])
    }
    # Read the file content
    lines <- readLines(module.file)
    # GBM is missing /// at the end
    if( lines[length(lines)] != "///" ){
        lines <- append(lines, "///")
    }
    
    keys <- c()
    modules <- list()
    module <- c()
    
    for (line in lines){
      
        if( line == "///" ){
        
            modules <- append(modules, list(module))
            module <- c()
        
        }else if( substr(line, 1, 1) == "M" ){
        
            key <- paste0(gsub("\t", " (", line), ")")
            keys <- append(keys, key)
        
        }else{
            module <- append(module, line)
        }
    }
    
    for( i in seq_along(modules) ){
        module.comp <- strsplit(modules[[i]], "\t")
        modules[[i]] <- lapply(module.comp, strsplit, split = ",")
    }
    
    names(modules) <- keys
    return(modules)
}
