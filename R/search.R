
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
    vertices <- lapply(
        sp$vpaths,
        function(path) V(graph)$name[path]
    )
    edges <- lapply(
        sp$epaths,
        function(path) E(graph)$source[path]
    )
    
    paths <- mapply(
        function(nodes, labels) {
            edges_str <- paste0(" -(", labels, ")-> ", nodes[-1])
            paste0(nodes[1], paste0(edges_str, collapse = ""))
        },
        vertices, edges
    )

    return(paths)
}
