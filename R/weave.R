
#' @name weavePath
#' @rdname weavePath

#' @importFrom igraph E<- k_shortest_paths as_data_frame
#' @importFrom MultiFactor MultiFactor
#' @importFrom KEGGREST keggLink
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

    if( to == "gbm" ){
        to <- from
        from <- "gbm"
    }
    
    FUN <- ifelse(from == "gbm", .path2gbm, .path2any)
    
    graph <- FUN(graph, from, to, k)
    # Add edge attribute to download critical edges
    E(graph)$mark <- E(graph)$mark != 0
    # Convert graph to data.frame
    graph_df <- as_data_frame(graph, what = "edges")
    # Select only edges to download
    graph_df <- graph_df[graph_df$mark, ]
    # Remove duplicate file paths
    graph_df <- graph_df[!duplicated(graph_df$path, incomparables = NA), ]
    # Initialise linkmap list
    linkmaps <- list()
    
    for( i in seq_len(nrow(graph_df)) ){

        if( graph_df$source[i] %in% c("ChocoPhlAn", "GBM", "GMM", "GO",
                                      "TIGRFAMs", "WoL")){
              
            cached <- .cache_resource(graph_df$path[i], graph_df$source[i])
            df <- read.csv(cached, colClasses = "character")
                
        }else if( graph_df$source[i] == "KEGG" ){
                
            kegg.link <- keggLink(graph_df$from[i], graph_df$to[i])
                
            df <- data.frame(
                x = gsub("^[^:]*:", "", kegg.link),
                y = gsub("^[^:]*:", "", names(kegg.link)),
                row.names = NULL
            )
        
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
