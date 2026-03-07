
#' @name searchPath
#' @rdname searchPath

#' @importFrom igraph E V k_shortest_paths
S7::method(searchPath, igraph) <- function(graph, by, k = 1){
    
    by.vars <- all.vars(by)
    
    from <- by.vars[1]
    to <- by.vars[2]
    
    # Find shortest path 
    sp <- k_shortest_paths(graph, from = from, to = to, k = k, mode = "all")
    # Convert vertex IDs to names
    vertices <- lapply(sp$vpaths, function(path) V(graph)$name[path])
    edges <- lapply(sp$epaths, function(path) E(graph)$source[path])
    
    paths <- mapply(
        function(nodes, labels) {
            edges_str <- paste0(" -(", labels, ")-> ", nodes[-1])
            paste0(nodes[1], paste0(edges_str, collapse = ""))
        },
        vertices, edges
    )

    return(paths)
}


#' @importFrom igraph k_shortest_paths E<- V<-
.draw_path <- function(graph, from, to, k){
    
    comb_df <- expand.grid(from = from, to = to, stringsAsFactors = FALSE)
    path.comb <- .generate_path_comb(k, nrow(comb_df))
    
    path_dfs <- list()
    
    for( i in seq_along(path.comb) ){
        
        j <- path.comb[[i]]
        orig <- comb_df$from[[i]]
        target <- comb_df$to[[i]]
        
        sp <- k_shortest_paths(
            graph, from = orig, to = target, k = j, mode = "all"
        )
        
        edge_idx <- sp$epaths[[j]]
        node_idx <- sp$vpaths[[j]]
        
        edges <- E(graph)$source[edge_idx]
        nodes <- V(graph)$name[node_idx]
        
        path_dfs[[i]] <- data.frame(
            from = nodes[-length(nodes)],
            to = nodes[-1],
            source = edges,
            step = seq(edges)
        )
    }
    # Bind paths
    path_df <- do.call(rbind, path_dfs)
    # Order by step
    path_df <- path_df[order(path_df$step), ]
    # Remove duplicates
    is.duplicate <- duplicated(path_df[ , c("from", "to", "source")])
    path_df <- path_df[!is.duplicate, ]
    return(path_df)
}
