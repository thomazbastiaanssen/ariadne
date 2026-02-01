
plot <- function(graph, by = NULL, k = 1){
  
    E(graph)$highlight <- FALSE
    E(graph)$name <- ""
    
    if( !is.null(by) ){
      
        by.vars <- all.vars(by)
        
        from <- by.vars[1]
        to <- by.vars[2]
        
        # Find shortest path 
        sp <- k_shortest_paths(graph, from = from, to = to, k = k, mode = "all")
        # Get edges IDs in path
        path_edges <- sp$epath[[k]]
        
        # Add edge attribute to highlight edges in the path
        E(graph)$highlight[path_edges] <- TRUE
        E(graph)$name[path_edges] <- E(graph)$source[path_edges]
    }

    # Plot graph with edges highlighted and others faded
    p <- ggraph(graph, layout = "stress") +
        geom_edge_link(
            aes(colour = highlight, label = name),
            edge_width = 1.2, fontface = "bold") +
        geom_node_point(size = 5, color = "darkorange") +
        geom_node_text(aes(label = name), vjust = 1.8, size = 4) +
        scale_edge_colour_manual(
            values = c(`TRUE` = "red", `FALSE` = "grey"),
            guide = "none"
        ) +
        theme_void()
    
    return(p)
}
