#' Draw path steps as a reproducible table
#' 
#' @name drawPath
#' 
#' @description
#' \code{drawPath} creates a data.frame with the metadata on each step in the
#' selected path, including its origin, version and url. This can be used to
#' complement analysis results and ensure reproducibility.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to plot. (Default: \code{NULL})
#' 
#' @param k \code{Numeric scalar}. The kth shortest path to plot.
#'   (Default: \code{1})
#' 
#' @param include \code{Character vector}. Nodes to cross in the path.
#'   (Default: \code{NULL})
#' 
#' @param exclude \code{Character vector}. Nodes to avoid in the path.
#'   (Default: \code{NULL})
#' 
#' @returns
#' A data.frame with path steps in the rows and step metadata in the columns.
#'
#' @examples
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Draw fourth path from ko to ec
#' df <- drawPath(graph, ko ~ ec, k = 4)
#' 
#' # View path metadata
#' df
NULL


#' @export
#' @rdname drawPath
setMethod("drawPath", signature = c(graph = "igraph"),
    function(graph, by, k = 1, include = NULL, exclude = NULL){
    # Draw path through graph
    path_df <- .draw_path(graph, by, k, include, exclude)
    # Add edge metadata for external use
    path_df <- .add_edge_metadata(path_df, graph, internal = FALSE)
    return(path_df)
})


.add_edge_metadata <- function(path_df, graph, internal){
    
    graph_df <- as_data_frame(graph, what = "both")
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
    
    versions <- attr(graph, "versions")
    path_df$version <- versions[path_df$source]
    
    graph_keys <- .get_edge_keys(edge_df)
    path_keys <- .get_edge_keys(path_df)
    
    idx <- match(path_keys, graph_keys)
    path_df$url <- edge_df$url[idx]
    
    if( internal ){
      
        path_df$initFrom <- path_df$from
        
        idy <- which(!is.na(path_df$url))
        path_df[idy, c("from", "to")] <- edge_df[idx[idy], c("from", "to")]
        
        # Retrieve specific names for KEGG, OTT and SPARQL queries
        path_df$specFrom <- .generic2specific(path_df, node_df, "from")
        path_df$specTo <- .generic2specific(path_df, node_df, "to")
    }
    
    return(path_df)
}
