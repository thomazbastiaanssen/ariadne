
#' @name weavePath
#' @rdname weavePath

#' @importFrom igraph E<- k_shortest_paths as_data_frame
#' @importFrom MultiFactor MultiFactor
#' @importFrom BiocParallel bplapply
S7::method(weavePath, igraph) <- function(graph, by, k = 1, init = NULL,
    prune = TRUE, verbose = TRUE, timeout = 1e6, ...){
    # Set timeout for downloads
    options(timeout = timeout)
    # Extract vars from formula
    by.vars <- all.vars(by)
    # Assign from and to vars
    from <- by.vars[1]
    to <- by.vars[2]
    
    graph_df <- as_data_frame(graph, what = "both")
    
    linkmaps <- list()
    path_df <- .draw_path(graph, from, to, k)
    
    for( i in seq_len(nrow(path_df)) ){
        # Retrieve step
        g <- path_df[i, ]
        # Print step
        if( verbose ) message(g$from, " -(", g$source, ")-> ", g$to)
        # Fetch linkmap
        linkmap <- .fetch_resource(
            graph_df, g$from, g$to, g$source, init, timeout, ...
        )
        # Add to linkmaps
        linkmaps[[paste0(g$from, "2", g$to)]] <- linkmap
        # Update init
        init <- if( prune ) unique(linkmap[[g$to]]) else NULL
    }
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    # Weave desired linkmap from MultiFactor
    linkmap <- weave(mf, by)
    return(linkmap)
}


#' @importFrom KEGGREST keggLink
#' @importFrom arrow read_parquet open_dataset
#' @importFrom dplyr filter collect
#' @importFrom stats na.omit
#' @importFrom rlang sym
.fetch_resource <- function(graph_df, from, to, repo, init, timeout, ...){
    
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
    
    is.source <- edge_df$source == repo
    
    idx <- which(
        edge_df$from == from & edge_df$to == to & is.source
    )
    
    if( length(idx) == 0L ){
        
        idx <- which(
            edge_df$from == to & edge_df$to == from & is.source
        )
    }
    
    if( length(idx) > 1L ){
        stop("Multiple matches were found.", call. = FALSE)
    }
    
    g <- edge_df[idx, , drop = FALSE]
    
    if( g$source %in% c("ChocoPhlAn", "GO", "TIGRFAMs", "WoL")){
        
        cached <- .cache_resource(g$path, g$source, c(g$from, g$to))
        
        if( !is.null(init) ){
            
            df <- cached |>
                open_dataset() |>
                filter(!!sym(from) %in% init) |>
                collect() |>
                as.data.frame()
        }else{
            df <- read_parquet(cached)
        }
        
    }else if( g$source == "KEGG" ){
    
        kegg.link <- keggLink(g$from, g$to)
        
        df <- data.frame(
            x = gsub("^[^:]*:", "", kegg.link),
            y = gsub("^[^:]*:", "", names(kegg.link)),
            row.names = NULL
        )
        
        if( !is.null(init) ){
            df <- df[df[[g$from]] %in% init, ]
        }
    
    }else if( g$source == "OTT" ){
        
        spec.from <- node_df[node_df$name == g$from, g$source]
        spec.to <- node_df[node_df$name == g$to, g$source]
        
        df <- data.frame(
            x = init,
            y = .queryOTT(init, spec.from, spec.to),
            row.names = NULL
        )
        
        df <- na.omit(df)
    
    }else if( g$source == "UniProt" ){
        
        g$from <- from
        g$to <- to
        
        spec.from <- node_df[node_df$name == g$from, g$source]
        spec.to <- node_df[node_df$name == g$to, g$source]
        
        spec.from <- ifelse(startsWith(spec.from, "uniref"), "uniref", spec.from)
        spec.to <- ifelse(startsWith(spec.to, "uniref"), "uniref", spec.to)
        
        df <- .querySPARQL(spec.from, spec.to, g$source, init, timeout, ...)
        df[] <- lapply(df, function(col) gsub("http.+/", "", col))
        
        if( spec.to == "uniref" ){
            df <- df[grepl(g$to, df[[spec.to]], ignore.case = TRUE), ]
            rownames(df) <- NULL
        }
    }
    # Add edge names
    colnames(df) <- c(g$from, g$to)
    return(df)
}
