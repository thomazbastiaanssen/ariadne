
#' @name searchPath
#' @rdname searchPath

#' @importFrom igraph E V k_shortest_paths
S7::method(searchPath, igraph) <- function(graph, by, k = 1){
    
    msg <- c()
    
    for( i in seq_len(k) ){
        
        msg <- c(msg, "Path ", i, ":\n")
        
        path_df <- .draw_path(graph, by, i, dup.rm = FALSE)
        
        unique_paths <- unique(path_df$path)
        
        for( j in unique_paths ){
            
            subpath_df <- path_df[path_df$path == j, ]
            
            path.str <- paste0(
                " -(", subpath_df$source, ")-> ", subpath_df$to, collapse = ""
            )
            
            path.str <- paste0(subpath_df$from[1], path.str)
            
            msg <- c(msg, path.str, "\n")
        }
        msg <- c(msg, "\n")
    }
    
    message(msg)
    invisible(NULL)
}


#' @importFrom stats ave
#' @importFrom igraph k_shortest_paths E<- V<- graph_from_data_frame topo_sort
.draw_path <- function(graph, by, k, dup.rm = TRUE){
    
    by.vars <- .formula2list(by)
    var.len <- length(by.vars)
    
    comb_df <- expand.grid(by.vars, stringsAsFactors = FALSE)
    path.comb <- .generate_path_comb(k, nrow(comb_df))
    
    path_dfs <- list()
        
    for( i in seq_len(nrow(comb_df)) ){
        
        orig <- comb_df[i, 1]
        target <- comb_df[i, ncol(comb_df)]
        cur.k <- path.comb[[i]]
        
        for( j in seq_len(ncol(comb_df) - 1) ){
            
            step.start <- comb_df[i, j]
            step.end <- comb_df[i, j + 1]
            
            sp <- k_shortest_paths(
                graph, step.start, step.end, k = cur.k, mode = "all"
            )
            
            edge_idx <- sp$epaths[[cur.k]]
            node_idx <- sp$vpaths[[cur.k]]
            
            edges <- E(graph)$source[edge_idx]
            nodes <- V(graph)$name[node_idx]
            
            path_dfs[[paste0(i, j)]] <- data.frame(
                from = nodes[-length(nodes)],
                to = nodes[-1],
                source = edges,
                path = paste0(orig, "2", target)
            )
        }
    }
    # Bind paths
    path_df <- do.call(rbind, path_dfs)
    # Add step number
    path_df$step <- ave(seq_along(path_df$path), path_df$path, FUN = seq_along)
    # If duplicates to be removed
    if( dup.rm ){
        # Remove duplicates
        is.duplicate <- duplicated(path_df[ , c("from", "to", "source")])
        path_df <- path_df[!is.duplicate, ]
        rownames(path_df) <- NULL
    }
    # Sort topology
    flow <- graph_from_data_frame(path_df)
    topo_order <- topo_sort(flow, mode = "out")
    # Get topological order
    node_order <- seq_along(topo_order)
    names(node_order) <- names(topo_order)
    # Order by graph topology
    node_order <- order(node_order[path_df$to])
    path_df <- path_df[node_order, ]
    return(path_df)
}


.expand_multiplex <- function(graph, by.vars){
  
    mod.name <- "gbm"
    
    mod.idx <- which(mod.name %in% by.vars)
    is.mod <- length(mod.idx) != 0L
    
    edge_df <- as_data_frame(graph, what = "edges")
    feat.name <- list(edge_df$to[edge_df$from == mod.name])
    
    if( is.mod ){
        after <- ifelse(mod.idx == length(by.vars), mod.idx - 1, mod.idx)
        by.vars <- append(by.vars, feat.name, after = after)
    }
    
    return(by.vars)
}


.generate_path_comb <- function(k, j){
    # Find minimum number of paths per neighbour
    num.paths <- rep(ceiling(k^(1 / j)), j)
    # Generate all combinations of path indices
    comb.indices <- expand.grid(lapply(num.paths, function(n) seq_len(n)))
    # Calculate the max coordinate per row
    comb.indices$max_val <- apply(comb.indices, 1L, max)
    
    var.cols <- setdiff(names(comb.indices), "max_val")
    comb.order <- do.call(order, as.list(comb.indices[, c("max_val", var.cols)]))
    # Order by max_val, then lex order
    comb.indices <- comb.indices[comb.order, ]
    rownames(comb.indices) <- NULL
    # Select the k-th combination of paths
    path.comb <- comb.indices[k, var.cols]
    return(path.comb)
}



.formula2list <- function(expr){
    parts <- .split_by_tilde(expr)
    lapply(parts, .split_by_plus)
}


.split_by_tilde <- function(expr) {
    if(is.call(expr) && expr[[1]] == as.name("~") ){
        # Recursively split left and right sides
        c(.split_by_tilde(expr[[2]]), .split_by_tilde(expr[[3]]))
    }else{
        # Base case: return expression as character
        list(expr)
    }
}


# Split each part by +
.split_by_plus <- function(expr){
    if( is.call(expr) && expr[[1]] == as.name("+") ){
        c(.split_by_plus(expr[[2]]), .split_by_plus(expr[[3]]))
    }else{
        trimws(deparse(expr))
    }
}

