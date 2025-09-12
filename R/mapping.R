

#' @rdname mapModules
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom anansi MultiFactor
mapModules <- function(mod.x, x.uniref, batch.size = 1000, drop.unmatched = FALSE){
    # Check arguments
    if( !(is.data.frame(mod.x) && ncol(mod.x) == 2) ){
        stop("'mod.x' must be a linkMap.", call. = FALSE)
    }
    if( !(is.data.frame(x.uniref) && ncol(x.uniref) == 2) ){
        stop("'x.uniref' must be a linkMap.", call. = FALSE)
    }
    if( !any(names(mod.x) %in% names(x.uniref)) ){
        stop("'mod.x' and 'x.uniref' should have one column in common.",
            call. = FALSE)
    }
    if( !any(names(x.uniref) %in% c("uniref50", "uniref90")) ){
        stop("One column of 'x.uniref' should contain uniref ids and be named ",
            "either uniref50 or uniref90.", call. = FALSE)
    }
    if( !is.numeric(batch.size) ){
        stop("'batch.size' must be a number.", call. = FALSE)
    }
    if( !is.logical(drop.unmatched) ){
        stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
    }
    # Keep only relevant bindings
    mf <- MultiFactor(list(a = mod.x, b = x.uniref), drop.unmatched = TRUE)
    # Retrieve unique uniref ids
    unique.uniref <- mf@levels[[which(colnames(mf) %in% c("uniref50", "uniref90"))]]
    # Split into query batches
    batches <- split(
        unique.uniref,
        ceiling(seq_along(unique.uniref) / batch.size)
    )
    # Query taxonomy from UniProt
    uniref.tax <- bplapply(batches, .querySPARQL)
    # Store taxa in modules list
    mf <- MultiFactor(
        list(a = mod.ko, b = ko.uniref90, c = uniref.tax),
        drop.unmatched = drop.unmatched
    )
    # Retrieve unique uniref ids
    
    return(mf)
}
