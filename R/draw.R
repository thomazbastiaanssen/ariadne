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
#' @param res.name \code{Character vector}. Names of resources to include in
#'   the graph. (Default: \code{NULL})
#' 
#' @param ... Unused.
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
    function(graph, by, k = 1, include = NULL, exclude = NULL, res.name = NULL){
    # Draw path through graph
    path_df <- .draw_path(graph, by, k, include, exclude, res.name)
    # Add edge metadata for external use
    path_df <- .add_edge_metadata(path_df, graph, internal = FALSE)
    return(path_df)
})


#' @importFrom igraph k_shortest_paths subgraph_from_edges E<- V<-
.draw_path <- function(graph, by, k, include, exclude, res.name,
    buffer.factor = 2, max.attempts = 5){
    # Check args
    if( !is.numeric(k) || length(k) != 1L || k <= 0 ){
        stop("'k' must be a positive integer.", call. = FALSE)
    }
    if( length(intersect(include, exclude)) != 0L ){
        stop("'include' and 'exclude' cannot overlap.", call. = FALSE)
    }
    # Extract by vars
    by.vars <- all.vars(by)
    from <- by.vars[1]
    to <- by.vars[2]
    # Subset graph by resource
    if( !is.null(res.name) ){
        graph <- subgraph_from_edges(graph, E(graph)[source %in% res.name])
    }
    # Initialise while vars
    j <- k
    i <- 0
    keep <- logical(0L)
    # Until enough paths found
    while( i < max.attempts && k > sum(keep) ){
        # Find paths
        sp <- k_shortest_paths(graph, from, to, k = j, mode = "all")
        # Select suitable paths
        keep <- vapply(
          sp$vpaths,
          function(v) all(include %in% names(v)) & !any(exclude %in% names(v)),
          logical(1L)
        )
        # Increase buffer and attempt
        j <- buffer.factor * j
        i <- i + 1
    }
    # Check results
    if( !any(keep) ){
        stop("No paths meet 'include' and 'exclude' criteria.", call. = FALSE)
    }
    if( k > sum(keep) ){
        stop("'k' is greater than the number of possible paths.", call. = FALSE)
    }
    # Select suitable paths
    sp <- lapply(sp, `[`, keep)
    # Find edges and nodes indices
    edge_idx <- sp$epaths[[k]]
    node_idx <- sp$vpaths[[k]]
    # Retrieve edges and nodes
    edges <- E(graph)$source[edge_idx]
    nodes <- V(graph)$name[node_idx]
    # Create path data.frame
    path_df <- data.frame(
        from = nodes[-length(nodes)],
        to = nodes[-1],
        source = edges
    )
    return(path_df)
}


#' @importFrom igraph as_data_frame
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
        # Store path variable order
        path_df$initFrom <- path_df$from
        path_df$initTo <- path_df$to
        # Use graph variable order for URL resources
        idy <- which(!is.na(path_df$url))
        path_df$from[idy] <- edge_df$from[idx[idy]]
        path_df$to[idy] <- edge_df$to[idx[idy]]
        # Retrieve specific names for the queries
        path_df$specFrom <- .generic2specific(path_df, node_df, "from")
        path_df$specTo <- .generic2specific(path_df, node_df, "to")
        path_df$specInitFrom <- .generic2specific(path_df, node_df, "initFrom")
    }
    return(path_df)
}
