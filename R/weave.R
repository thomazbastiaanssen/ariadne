
#' @name weavePath
#' @rdname weavePath

#' @importFrom igraph E<- k_shortest_paths as_data_frame
#' @importFrom MultiFactor MultiFactor
#' @importFrom BiocParallel bplapply
S7::method(weavePath, igraph) <- function(graph, by, k = 1, timeout = 1e6){
    # Set timeout for downloads
    options(timeout = timeout)
    # Extract vars from formula
    by.vars <- all.vars(by)
    # Assign from and to vars
    from <- by.vars[1]
    to <- by.vars[2]
    
    # FUN <- ifelse("gbm" %in% by.vars, .path2gbm, .path2any)
    # graph <- FUN(graph, from, to, k)
    
    graph_df <- as_data_frame(graph, what = "both")
    
    linkmaps <- list()
    path_df <- .draw_path(graph, from, to, k)
    
    for( i in seq_len(nrow(path_df)) ){
      # Define step vars
      x <- path_df[i, "from"]
      y <- path_df[i, "to"]
      z <- path_df[i, "source"]
      # Fetch linkmap
      linkmap <- .fetch_resource(graph_df, x, y, z) # Add prev linkmap as input
      # Add to linkmaps
      linkmaps[[paste0(x, "2", y)]] <- linkmap
    }
    
    
    return(linkmaps)
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    # Weave desired linkmap from MultiFactor
    linkmap <- weave(mf, by)
    return(linkmap)
}



#' @importFrom KEGGREST keggLink
#' @importFrom arrow read_parquet
.fetch_resource <- function(graph_df, from, to, repo){
    
    edge_df <- graph_df$edges
    
    idx <- which(
        edge_df$from == from & edge_df$to == to & edge_df$source == repo
    )
    
    if( length(idx) == 0L ){
        tmp <- to
        to <- from
        from <- tmp
        
        idx <- which(
            edge_df$from == from & edge_df$to == to & edge_df$source == repo
        )
    }
    
    if( length(idx) > 1L ){
        stop("Multiple matches were found.", call. = FALSE)
    }
    
    g <- edge_df[idx, ]
    
    if( g$source %in% c("ChocoPhlAn", "GM", "GO", "TIGRFAMs", "WoL")){
    
        cached <- .cache_resource(g$path, g$source)
        df <- read_parquet(cached)
    
    }else if( g$source == "KEGG" ){
    
        kegg.link <- keggLink(g$from, g$to)
        
        df <- data.frame(
            x = gsub("^[^:]*:", "", kegg.link),
            y = gsub("^[^:]*:", "", names(kegg.link)),
            row.names = NULL
        )
    
    }else if( g$source == "UniProt" ){
        # Proceed with UniProt queries
        linkmap <- .querySPARQL(x, g$from, g$to, g$source[i])
        linkmap <- apply(
            linkmap, 2L, function(x) gsub("([^/]+)$", "", x), simplify = FALSE
        )
    }
    
    # Add edge names
    colnames(df) <- c(g$from, g$to)
    return(df)
}
