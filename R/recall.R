#' Reproduce path from a path data.frame
#' 
#' @name recallPath
#' @aliases recallPath
#' 
#' @description
#' \code{recallPath}
#' 
#' @param path_df A data.frame object.
#' 
#' @param mode \code{Character scalar}. A string specifying the weaving mode,
#'   either \code{"simple"} or \code{"complex"}. (Default: \code{"simple"})
#' 
#' @param ... Additional arguments.
#' 
#' @returns
#' A two-column data.frame (x-to-y linkmap), with optional names as a third
#' column.
#' 
#' @examples
#' # Import example path df
#' data("pathMeta", package = "ariadne")
#' 
#' # Weave path starting from initial values
#' tax2bugsig <- recallPath(pathMeta, init = tax.labs)
#' 
#' head(tax2bugsig)
NULL


#' @export
#' @rdname recallPath
#' @importFrom stats as.formula
#' @importFrom igraph as_data_frame subgraph_from_edges
setMethod("recallPath", signature = c(path_df = "data.frame"),
    function(path_df, mode = "simple", ...){
    
    if( !mode %in% c("simple", "complex") ){
        stop("'mode' must be either 'simple' or 'complex'.", call. = FALSE)
    }
    
    versions <- as.list(path_df$version)
    names(versions) <- path_df$source
    
    graph <- ariadne(versions = versions)
    
    E(graph)$name <- graph |>
        as_data_frame(what = "edges") |>
        .get_edge_keys()
    
    keep <- .get_edge_keys(path_df)
    
    graph <- subgraph_from_edges(graph, keep)
    
    by <- c(path_df$from[1L], path_df$to[nrow(path_df)]) |>
        paste(collapse = "~") |>
        as.formula()
    
    FUN <- switch(mode, simple = weavePath, complex = weaveComplex)
    
    linkmap <- FUN(graph, by, ...)
    return(linkmap)
})
