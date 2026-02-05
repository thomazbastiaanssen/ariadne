
#' @name plotPath
#' @rdname plotPath

#' @importFrom igraph E k_shortest_paths
#' @importFrom ggplot2 aes theme_void
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_colour_manual
S7::method(plotPath, igraph) <- function(graph, by = NULL, k = 1){
    
    E(graph)$mark <- 0
    E(graph)$name <- ""
    
    if( !is.null(by) ){
      
        by.vars <- all.vars(by)
        
        from <- by.vars[1]
        to <- by.vars[2]
        
        if( to == "gbm" ){
            to <- from
            from <- "gbm"
        }
        
        FUN <- ifelse(from == "gbm", .path2gbm, .path2any)
        
        graph <- FUN(graph, from, to, k)
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
    # Include grey for edges not in paths
    path.colours <- c(path.colours, "0" = "grey80")
    # Plot graph with edges marked and others faded
    p <- ggraph(graph, layout = "stress") +
        geom_edge_link(
            aes(colour = factor(mark), label = name),
            edge_width = 1.2, fontface = "bold",
            show.legend = TRUE) +
        geom_node_point(size = 5, color = "darkorange") +
        geom_node_text(aes(label = name), vjust = 1.8, size = 4) +
        scale_edge_color_manual(
            values = path.colours,
            breaks = path.names,
            labels = path.names,
            name = "Paths"
        ) +
        theme_void() +
        theme(legend.position = "bottom")
    
    return(p)
}

.path2any <- function(graph, from, to, k){
    # Find shortest path 
    sp <- k_shortest_paths(graph, from = from, to = to, k = k, mode = "all")
    # Get edges IDs in path
    path.edges <- unique(unlist(sp$epaths[1:k]))
    # Assign positive mark to path edges
    E(graph)$mark[path.edges] <- 1
    # Return graph with aesthetics defined
    return(graph)
}


.path2gbm <- function(graph, from, to, k) {
    # Find unique neighbours
    neighbours <- unique(neighbors(graph, from, mode = "all"))
    # Gather all paths for each intermediate node up to k
    paths <- list()
    for (i in seq_along(neighbours)) {
        sp <- k_shortest_paths(
            graph, from = neighbours[i], to = to, k = k, mode = "all"
        )
        paths[[i]] <- sp$epaths[1:min(k, length(sp$epaths))]
    }
    # Number of paths per neighbor
    num.paths <- lengths(paths)
    # Generate all combinations of path indices
    comb.indices <- expand.grid(lapply(num.paths, function(n) seq_len(n)))
    # Calculate the max coordinate per row
    comb.indices$max_val <- apply(comb.indices, 1, max)
    
    var_cols <- setdiff(names(comb.indices), "max_val")
    comb.order <- do.call(order, as.list(comb.indices[, c("max_val", var_cols)]))
    
    # Order by max_val, then lex order
    comb.indices <- comb.indices[comb.order, ]
    # Select the k-th combination of paths
    chosen.comb <- comb.indices[k, var_cols]
    # Assign unique mark per neighbour to colour edges accordingly
    for (i in seq_along(chosen.comb)) {
        path.edges <- unlist(paths[[i]][[chosen.comb[[i]]]])
        # Assign mark to these edges
        E(graph)$mark[path.edges] <- i
    }
    # Add edges from 'from' to intermediate nodes (original edges)
    from_id <- which(V(graph)$name == from)
    neighbour.ids <- as.integer(neighbours)
    
    for (i in seq_along(neighbours)) {
        eid <- get_edge_ids(graph, c(from_id, neighbour.ids[i]))
        if (eid != 0) {
            E(graph)$mark[eid] <- i
        }
    }
    
    return(graph)
}

