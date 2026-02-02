
#' @name ariadne
#' @rdname ariadne

S7::method(ariadne, S7::class_list) <-
    function(resources, versions = NULL){
    # Convert resources to graph data.frames
    graph_dfs <- mapply(
        .resource2graph_df,
        resources, names(resources),
        MoreArgs = list(versions = versions),
        SIMPLIFY = FALSE
    )
    # Merge graph data.frames 
    graph_df <- do.call(rbind, graph_dfs)
    rownames(graph_df) <- NULL
    # Create igraph object from graph data.frame
    graph <- graph_from_data_frame(graph_df, directed = TRUE)
    return(graph)
}


.resource2graph_df <- function(res, name, versions){
    # Get existing feature combinations
    graph_df <- expand.grid(
        from = names(res$from), to = names(res$to), stringsAsFactors = FALSE
    )
    # Add source name
    graph_df$source <- name
    # Select resource version
    if(name %in% names(versions)){
        version <- versions[[name]]
    }else{
        version <- res$version
    }
    # Construct file path
    if( is.null(res$path) ){
        graph_df$path <- NA
    } else {
        graph_df$path <- res$path(res$repo, version, graph_df$from, graph_df$to)
    }
    return(graph_df)
}
