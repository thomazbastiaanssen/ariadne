#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr bind_rows
mapModules <- function(modules, map){
    # Keep only relevant bindings
    keep <- names(map) %in% unique(unlist(modules))
    map <- map[keep]
    # Query taxonomy from UniRef90
    tax.list <- bplapply(map, .querySPARQL)
    # Store taxa in modules list
    sig.list <- vector("list", length(modules))
    for( i in seq_along(modules) ){
        # Keep only relevant bindings
        keep <- names(tax.list) %in% modules[[i]]
        sig.list[[i]] <- unname(unlist(tax.list[keep]))
    }
    # Add module names
    names(sig.list) <- names(modules)
    return(sig.list)
}
