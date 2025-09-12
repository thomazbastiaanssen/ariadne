#' @export
#' @importFrom BiocParallel bplapply
mapModules <- function(modules, map, remove.empty = TRUE, message = TRUE){
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
    keep <- names(map) %in% unique(unlist(modules, use.names = FALSE))
    map <- map[keep]
    if( message ){
        message("Mapping ", length(modules), " modules to ",
            length(unlist(map)), " uniref ids.")
    }
    # Query taxonomy from UniRef90
    tax.list <- bplapply(map, .querySPARQL)
    # Find functions for each taxon
    tax.linkmap <- as.linkMap(tax.list)
    tax <- split(tax.linkmap$x, tax.linkmap$y)
    if( message ){
        message(length(tax), " tax ids were retrieved from UniProt.")
    }
    # Store taxa in modules list
    sig.list <- bplapply(modules, function(module) {
        members <- vapply(tax, function(tax.item) {
            all(vapply(module, function(mod) any(mod %in% tax.item), logical(1)))
        }, logical(1))
        names(tax)[members]
    })
    if( remove.empty ){
        # Remove empty modules 
        sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
    }
    return(sig.list)
}

as.linkMap <- function(values, keys = NULL, col.names = NULL){
    if( is.null(keys) ){
        keys <- names(values)
    }
    # Create linkMap
    linkmap <- data.frame(
        x = rep(keys, vapply(values, length, 1)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    # Assign custom colnames
    if( !is.null(col.names) ){
        names(linkmap) <- col.names
    }
    return(linkmap)
}
