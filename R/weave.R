
#' @name weavePath
#' @rdname weavePath

#' @importFrom methods getClass
#' @importFrom igraph k_shortest_paths
S7::method(weavePath, S7::class_any) <-
    function(graph, by, k = 1, timeout = 1e6){
    
    # Set timeout for downloads
    options(timeout = timeout)
    # Extract vars from formula
    by.vars <- all.vars(by)
    # Assign from and to vars
    from <- by.vars[1]
    to <- by.vars[2]
    # Find shortest path from A to C
    sp <- k_shortest_paths(graph, from = from, to = to, k = k, mode = "all")
    path_edges <- sp$epaths[[k]]  # edge IDs in path
    # Add edge attribute to download critical edges
    E(graph)$download <- FALSE
    E(graph)$download[path_edges] <- TRUE
    # Convert graph to data.frame
    graph_df <- as_data_frame(graph, what = "edges")
    # Select only edges to download
    graph_df <- graph_df[graph_df$download, ]
    # Initialise linkmap list
    linkmaps <- list()
    
    for( i in seq_len(nrow(graph_df)) ){

        if( graph_df$source[i] %in% c("ChocoPhlAn", "WoL", "GO")){
              
            cached <- .cache_resource(graph_df$path[i], graph_df$source[i])
            df <- read.csv(cached)
                
        }else if( graph_df$source[i] == "KEGG" ){
                
            kegg.link <- keggLink(graph_df$from[i], graph_df$to[i])
                
            df <- data.frame(
                x = gsub("^[^:]*:", "", kegg.link),
                y = gsub("^[^:]*:", "", names(kegg.link)),
                row.names = NULL
            )
        
        }else if( graph_df$source[i] %in% c("GBM", "GMM") ){
            
            
        }else if( graph_df$source[i] == "UniProt" ){
              
            # .querySPARQL(with uniref features from previous step)
              
        }
        # Add edge names
        colnames(df) <- c(graph_df$from[i], graph_df$to[i])
        # Append linkmap
        linkmaps[[length(linkmaps) + 1]] <- df
    }
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    # Weave desired linkmap from MultiFactor
    linkmap <- weave(mf, by)
    return(linkmap)
}
