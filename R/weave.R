
#' @name weavePath
#' @rdname weavePath

#' @importFrom igraph E<- k_shortest_paths as_data_frame
#' @importFrom MultiFactor MultiFactor
#' @importFrom BiocParallel bplapply
S7::method(weavePath, igraph) <- function(graph, by, k = 1, timeout = 1e6){
    # Set timeout for downloads
    options(timeout = timeout)
    # Initialise mark for files to download
    E(graph)$mark <- 0
    # Extract vars from formula
    by.vars <- all.vars(by)
    # Assign from and to vars
    from <- by.vars[1]
    to <- by.vars[2]
    
    FUN <- ifelse("gbm" %in% by.vars, .path2gbm, .path2any)
    graph <- FUN(graph, from, to, k)
    # Add edge attribute to download critical edges
    E(graph)$mark <- E(graph)$mark != 0
    # Convert graph to data.frame
    graph_df <- as_data_frame(graph, what = "edges")
    # Select only edges to download
    graph_df <- graph_df[graph_df$mark, , drop = FALSE]
    # Remove duplicate file paths
    graph_df <- graph_df[!duplicated(graph_df$path, incomparables = NA), ]
    # Split by UniProt which needs past linkmaps
    uniprot_idx <- graph_df$source == "UniProt"
    graph_df1 <- graph_df[!uniprot_idx, , drop = FALSE]
    graph_df2 <- graph_df[uniprot_idx, , drop = FALSE]
    # Fetch resources in parallel
    linkmaps <- lapply(
        seq_len(nrow(graph_df1)),
        function(i) .fetch_resource(graph_df1[i, ])
    )
    # Add names to linkmaps
    names(linkmaps) <- paste0(graph_df1$from, "2", graph_df1$to)
    # Proceed with UniProt queries
    for( i in seq_len(nrow(graph_df2)) ){
    
        input.ind <- grepl(paste0(".+2", graph_df2$from[i]), names(linkmaps))
        x <- unique(linkmaps[input.ind][[1]][ , graph_df2$from[i]])
        
        from <- ifelse(grepl("uniref", graph_df2$from[i]), "uniref", graph_df2$from[i])
        to <- ifelse(grepl("uniref", graph_df2$to[i]), "uniref", graph_df2$to[i])
        
        df <- .querySPARQL(x, from, to, graph_df2$source[i])
        df <- apply(df, 2L, function(x) gsub("([^/]+)$", "", x), simplify = FALSE)
        
        # Add edge names
        colnames(df) <- c(graph_df2$from[i], graph_df2$to[i])
        # Append linkmap
        linkmaps[[paste0(graph_df2$from[i], "2", graph_df2$to[i])]] <- df
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
.fetch_resource <- function(g){
    
    if( g$source %in% c("ChocoPhlAn", "GBM", "GMM", "GO", "TIGRFAMs", "WoL")){
    
        cached <- .cache_resource(g$path, g$source)
        df <- read_parquet(cached)
    
    }else if( g$source == "KEGG" ){
    
        kegg.link <- keggLink(g$from, g$to)
        
        df <- data.frame(
            x = gsub("^[^:]*:", "", kegg.link),
            y = gsub("^[^:]*:", "", names(kegg.link)),
            row.names = NULL
        )
    
    }
    
    # Add edge names
    colnames(df) <- c(g$from, g$to)
    return(df)
}
