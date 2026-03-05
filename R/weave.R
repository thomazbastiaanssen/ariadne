
#' @name weavePath
#' @rdname weavePath

#' @importFrom igraph E<- k_shortest_paths as_data_frame
#' @importFrom MultiFactor MultiFactor
#' @importFrom BiocParallel bplapply
S7::method(weavePath, igraph) <- function(
    graph, by, k = 1, init = NULL, prune = TRUE, verbose = TRUE, timeout = 1e6
    ){
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
        # Print step
        if( verbose ) message(x, " -(", z, ")-> ", y)
        # Fetch linkmap
        linkmap <- .fetch_resource(graph_df, x, y, z, init, timeout)
        # Add to linkmaps
        linkmaps[[paste0(x, "2", y)]] <- linkmap
        # Update init
        init <- if( prune ) linkmap[[y]] else NULL
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
.fetch_resource <- function(graph_df, from, to, repo, init, timeout){
    
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
    
    if( g$source %in% c("ChocoPhlAn", "GM", "GO", "TIGRFAMs", "WoL")){
        
        cached <- .cache_resource(g$path, g$source, c(g$from, g$to))
        
        if( !is.null(init) ){
            
            df <- cached |>
                open_dataset() |>
                filter(!!sym(from) %in% init) |>
                collect()
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
        
        df <- .querySPARQL(spec.from, spec.to, g$source, init, timeout)
        df[] <- lapply(df, function(col) gsub(".*/", "", col))
        
        if( spec.to == "uniref" ){
            df <- df[grepl(g$to, df[[spec.to]], ignore.case = TRUE), ]
        }
    }
    # Add edge names
    colnames(df) <- c(g$from, g$to)
    return(df)
}
