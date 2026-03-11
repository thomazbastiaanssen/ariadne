#' Plot paths between resources
#'
#' @name plotPath
#' @rdname plotPath
#' 
#' @description
#' \code{plotPath} provides a visual of a graph and the selected paths.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to plot. (Default: \code{NULL})
#' 
#' @param k \code{Numeric scalar}. The kth shortest path to plot.
#'   (Default: \code{1})
#' 
#' @param focus \code{Logical scalar}. Whether the edges and nodes in the path
#'   should be plotted. (Default: \code{FALSE})
#' 
#' @return
#' A ggplot2 object
#' 
#' @examples
#' 
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Plot fifth path from ko to ec
#' plotPath(graph, ko ~ ec, k = 5)
#' 
#' # Plot first path through uniref90
#' plotPath(graph, taxname ~ uniref90 ~ ko)
#' 
#' # Plot first path from ko and ec to uniref90
#' plotPath(graph, ko + ec ~ uniref90)
#' 
NULL

#' @importFrom igraph as_data_frame graph_from_data_frame ends E<- V<-
#' @importFrom ggplot2 aes scale_alpha theme_void theme
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_colour_manual
S7::method(plotPath, igraph) <- function(graph, by = NULL, k = 1, focus = FALSE){
    # Check args
    if( !is.numeric(k) || length(k) != 1L || k <= 0 ){
        stop("'k' must be a positive integer.", call. = FALSE)
    }
    if( focus && is.null(by) ){
        stop("'focus' can be TRUE when 'by' is defined.", call. = FALSE)
    }
    
    graph_df <- as_data_frame(graph, what = "both")
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
  
    edge_df$mark <- 0
    edge_df$name <- ""
    edge_df$alpha <- 1
    node_df$alpha <- 1
    
    path_names <- character(0L)
    alpha_min <- 1
    
    if( !is.null(by) ){
        
        path_df <- .draw_path(graph, by, k)
        
        graph_keys <- .get_edge_keys(edge_df)
        path_keys <- .get_edge_keys(path_df)
        
        path_names <- unique(path_df$path)
        
        edge_df$mark <- path_df$path[match(graph_keys, path_keys)]
        edge_df$mark[is.na(edge_df$mark)] <- 0
    }
    # Add edge attribute to mark edges in the path
    keep <- edge_df$mark != 0
    edge_df$name[keep] <- edge_df$source[keep]
    # Define path colours
    path_colours <- rainbow(length(path_names))
    names(path_colours) <- path_names
    # Include grey for edges not in paths
    path_colours <- c(path_colours, "0" = "grey80")
    # Create graph from edges and nodes data
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    # Create a vector for edge alpha: 1 if marked, else 0 (transparent)
    if( focus ){
        E(graph)$alpha <- ifelse(E(graph)$mark != 0, 1, 0)
        connected_nodes <- unique(c(ends(graph, E(graph)[E(graph)$mark != 0])))
        V(graph)$alpha <- ifelse(V(graph)$name %in% connected_nodes, 1, 0)
        alpha_min <- 0
    }
    # Plot graph with edges marked and others faded
    p <- ggraph(graph, layout = "stress") +
        geom_edge_link(
            aes(colour = factor(mark), label = name, alpha = alpha),
            edge_width = 1.2, fontface = "bold",
            show.legend = TRUE
        ) +
        geom_node_point(aes(alpha = alpha), size = 5, colour = "darkorange") +
        geom_node_text(aes(label = name, alpha = alpha), vjust = 1.8, size = 4) +
        scale_edge_colour_manual(
            values = path_colours,
            breaks = path_names,
            labels = path_names,
            name = "Paths"
        ) +
        scale_alpha(range = c(alpha_min, 1), guide = "none") +
        theme_void() +
        theme(legend.position = "bottom")
    
    return(p)
}


# Get edge keys from edges data
.get_edge_keys <- function(g) {
    edges <- mapply(
        function(x, y, z) paste(c(sort(c(x, y)), z), collapse = "_"),
        x = g$from, y = g$to, z = g$source,
        USE.NAMES = FALSE
    )
    return(edges)
}

