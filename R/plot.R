#' Plot paths between resources
#' 
#' @name plotPath
#' 
#' @description
#' \code{plotPath} provides a visual of a graph and the selected path.
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
#' @param focus \code{Logical scalar}. Whether the edges and nodes in the path
#'   should be plotted. (Default: \code{FALSE})
#' 
#' @param ... Unused.
#' 
#' @returns A ggplot2 object.
#' 
#' @examples
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Plot graph for a subset of resources
#' plotPath(graph, res.name = c("KEGG", "WoL"))
#' 
#' # Plot fifth path from ko to ec
#' plotPath(graph, ko ~ ec, k = 5)
#' 
#' # Plot first path including uniref90
#' plotPath(graph, taxname ~ ko, include = "uniref90")
#' 
#' # Plot first 5 paths excluding uniref50 and uniref100
#' plotPath(graph, taxname ~ ko, k = 5, exclude = "uniref50")
NULL

#' @export
#' @rdname plotPath
#' @importFrom igraph as_data_frame graph_from_data_frame ends
#' @importFrom ggplot2 aes theme_void theme
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_colour_manual
setMethod("plotPath", signature = c(graph = "igraph"),
    function(graph, by = NULL, k = 1, include = NULL, exclude = NULL,
    res.name = NULL, focus = FALSE){
    # Check args
    if( !is.logical(focus) || length(focus) != 1L ){
        stop("'focus' must be TRUE or FALSE.", call. = FALSE)
    }
    
    by_null <- is.null(by)
    
    if( focus && by_null ){
        stop("'focus' can be TRUE when 'by' is defined.", call. = FALSE)
    }
    
    graph_df <- as_data_frame(graph, what = "both")
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
  
    edge_df$mark <- 0
    edge_df$name <- ""
    edge_df$alpha <- TRUE
    node_df$alpha <- TRUE

    if( !by_null ){
        
        path_df <- .draw_path(graph, by, k, include, exclude, res.name)
        
        graph_keys <- .get_edge_keys(edge_df)
        path_keys <- .get_edge_keys(path_df)
        
        edge_df$mark <- as.integer(graph_keys %in% path_keys)
    }
    # Add edge attribute to mark edges in the path
    keep <- edge_df$mark != 0
    edge_df$name[keep] <- edge_df$source[keep]
    # Include grey for edges not in paths
    path_colours <- c("0" = "grey80", "1" = "red")
    # Create a vector for edge alpha: 1 if marked, else 0 (transparent)
    if( focus ){
        edge_df$alpha <- edge_df$mark != 0
        connected_nodes <- unique(c(ends(graph, E(graph)[edge_df$mark != 0])))
        node_df$alpha <- node_df$name %in% connected_nodes
    }else if( !is.null(res.name) ){
        edge_df$alpha <- edge_df$source %in% res.name
        node_df$alpha <- rowSums(!is.na(node_df[ , res.name, drop = FALSE])) != 0L
    }
    # Create graph from edges and nodes data
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    # Plot graph with edges marked and others faded
    p <- ggraph(graph, layout = "stress") +
        geom_edge_link(aes(colour = factor(.data$mark), label = .data$name,
            alpha = .data$alpha), edge_width = 1, fontface = "bold",
            show.legend = TRUE) +
        geom_node_point(aes(filter = .data$alpha),
            size = 4, colour = "darkorange") +
        geom_node_text(aes(label = .data$name, filter = .data$alpha),
            vjust = 1.8, size = 4) +
        scale_edge_colour_manual(values = path_colours) +
        theme_void() +
        theme(legend.position = "none")
    
    return(p)
})


# Get edge keys from edges data
.get_edge_keys <- function(g) {
    edges <- mapply(
        function(x, y, z) paste(c(sort(c(x, y)), z), collapse = "_"),
        x = g$from, y = g$to, z = g$source,
        USE.NAMES = FALSE
    )
    return(edges)
}

