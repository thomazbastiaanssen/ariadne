#' Search path between resources
#' 
#' @name searchPath
#' @rdname searchPath
#' 
#' @description
#' \code{searchPath} allows to search the first k-th shortest paths between
#' resources.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to search.
#' 
#' @param k \code{Numeric scalar}. The kth shortest paths to search.
#'   (Default: \code{1})
#' 
#' @param include \code{Character vector}. Nodes to cross in the path.
#'   (Default: \code{NULL})
#' 
#' @param exclude \code{Character vector}. Nodes to avoid in the path.
#'   (Default: \code{NULL})
#' 
#' @return
#' NULL, message
#' 
#' @examples
#' 
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Search first 5 paths from ko to ec
#' searchPath(graph, ko ~ ec, k = 5)
#' 
#' # Search first path including uniref90
#' searchPath(graph, taxname ~ ko, include = "uniref90")
#' 
#' # Search first 5 paths excluding uniref50 and uniref100
#' searchPath(graph, taxname ~ ko, k = 5, exclude = c("uniref50", "uniref100"))
#' 
NULL


#' @importFrom igraph E V k_shortest_paths
S7::method(searchPath, igraph) <- function(
    graph, by, k = 1, include = NULL, exclude = NULL){
    # Initialise message
    msg <- c()
    # Print paths up to k
    for( j in seq_len(k) ){
        # Add path number
        msg <- c(msg, "Path ", j, ":\n")
        # Get path
        path_df <- .draw_path(graph, by, j, include, exclude)
        # Add path string
        path_str <- paste0(
            " -(", path_df$source, ")-> ", path_df$to, collapse = ""
        )
        # Add origin
        path_str <- paste0(path_df$from[1], path_str)
        # Add new line 
        msg <- c(msg, path_str, "\n\n")
    }
    # Send message
    message(msg)
    invisible(NULL)
}


#' @importFrom igraph k_shortest_paths E<- V<-
.draw_path <- function(
    graph, by, k, include, exclude, buffer.factor = 2, max.attempts = 5){
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


.exclude_gm <- function(graph, by.vars){
    
    mod.names <- c("gbm", "gmm")
    # Identify module name
    is.mod <- by.vars %in% mod.names
    
    mod.idx <- which(mod.name %in% by.vars)
    is.mod <- length(mod.idx) != 0L
    
    edge_df <- as_data_frame(graph, what = "edges")
    feat.name <- paste(edge_df$to[edge_df$from == mod.name], collapse = "+")
    
    if( is.mod ){
        after <- ifelse(mod.idx == length(by.vars), mod.idx - 1, mod.idx)
        by.vars <- append(by.vars, feat.name, after = after)
    }
    
    return(by.vars)
}


.generate_path_comb <- function(k, j){
    # Find minimum number of paths per neighbour
    num_paths <- rep(ceiling(k^(1 / j)), j)
    # Generate all combinations of path indices
    comb_indices <- expand.grid(lapply(num_paths, function(n) seq_len(n)))
    # Calculate the max and total per row
    comb_points <- list(apply(comb_indices, 1L, max), rowSums(comb_indices))
    # Sort indices by max and sum values
    comb_indices <- comb_indices[do.call(order, comb_points), , drop = FALSE]
    rownames(comb_indices) <- NULL
    # Select the k-th combination of paths
    path_comb <- comb_indices[k, ]
    return(path_comb)
}

