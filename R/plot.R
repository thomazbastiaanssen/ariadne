
#' @name plotPath
#' @rdname plotPath

#' @importFrom igraph as_data_frame graph_from_data_frame ends E<- V<-
#' @importFrom ggplot2 aes scale_alpha theme_void theme
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_colour_manual
S7::method(plotPath, igraph) <- function(graph, by = NULL, k = 1, focus = FALSE){
    if( focus && is.null(by) ){
        stop("'focus' can be TRUE when 'by' is defined.", call. = FALSE)
    }
  
    E(graph)$mark <- 0
    E(graph)$name <- ""
    E(graph)$alpha <- 1
    V(graph)$alpha <- 1
    path_levels <- ""
    alpha_min <- 1
    
    if( !is.null(by) ){
        
        path_df <- .draw_path(graph, by, k)
        
        graph_keys <- .get_edge_keys(graph)
        path_keys <- .get_edge_keys(graph_from_data_frame(path_df))

        path_levels <- factor(path_df$path, levels = unique(path_df$path))
        
        key2level <- as.integer(path_levels)
        names(key2level) <- path_keys
        
        E(graph)$mark <- key2level[graph_keys]
        E(graph)$mark[is.na(E(graph)$mark)] <- 0
        
    }
    # Add edge attribute to mark edges in the path
    keep <- E(graph)$mark != 0
    E(graph)$name[keep] <- E(graph)$source[keep]
    # Define path names
    num.paths <- length(unique(E(graph)$mark))
    path.names <- as.character(seq_len(num.paths))
    # Define path colours
    path.colours <- rainbow(num.paths)
    names(path.colours) <- path.names
    # Create a vector for edge alpha: 1 if marked, else 0 (transparent)
    if( focus ){
        E(graph)$alpha <- ifelse(E(graph)$mark != 0, 1, 0)
        connected_nodes <- unique(c(ends(graph, E(graph)[E(graph)$mark != 0])))
        V(graph)$alpha <- ifelse(V(graph)$name %in% connected_nodes, 1, 0)
        alpha_min <- 0
    }
    # Include grey for edges not in paths
    path.colours <- c(path.colours, "0" = "grey80")
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
            values = path.colours,
            breaks = path.names,
            labels = c(levels(path_levels), ""),
            name = "Paths"
        ) +
        scale_alpha(range = c(alpha_min, 1), guide = "none") +
        theme_void() +
        theme(legend.position = "bottom")
    
    return(p)
}


# Function to get sorted edge keys for a graph
.get_edge_keys <- function(graph) {
    edges <- as_edgelist(graph)
    repo <- E(graph)$source
    edges <- apply(cbind(edges, repo), 1L, function(row) {
        vertices <- sort(row[1:2])
        paste(c(vertices, row[3]), collapse = "_")
    })
    return(edges)
}

