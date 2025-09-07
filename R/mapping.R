#' @export
#' @importFrom BiocParallel bplapply
mapModules <- function(modules, map, remove.empty = FALSE){
    # Check arguments
    if( !is.vector(modules) ){
      stop("'modules' must be a character vector or list of character ",
          "vectors, where each vector corresponds to a module.", call. = FALSE)
    }
    if( !is.vector(map) ){
      stop("'map' must be a character vector or list of character ",
          "vectors, where each vector corresponds to a mapping.", call. = FALSE)
    }
    if( !is.logical(remove.empty) ){
        stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
    }
    # Keep only relevant bindings
    keep <- names(map) %in% unique(unlist(modules))
    map <- map[keep]
    # Query taxonomy from UniRef90
    tax.list <- bplapply(map, .querySPARQL)
    # Store taxa in modules list
    sig.list <- bplapply(modules, function(module){
        keep <- names(tax.list) %in% module
        module <- unname(unlist(tax.list[keep]))
        return(module)
    })
    if( remove.empty ){
        # Remove empty modules 
        sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
    }
    return(sig.list)
}
