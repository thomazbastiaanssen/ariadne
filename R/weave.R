#' Weave path between resources
#' 
#' @name weavePath
#' @rdname weavePath
#' 
#' @description
#' \code{weavePath} weaves the path between resources.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to weave.
#' 
#' @param k \code{Numeric scalar}. The kth shortest path to weave.
#'   (Default: \code{1})
#' 
#' @param init \code{Character vector}. Initial values to prune the first
#'   mapping step. (Default: \code{NULL})
#' 
#' @param prune \code{Logical scalar}. Should the values from each step be used
#'   prune the following step. (Default: \code{TRUE})
#' 
#' @param output.format \code{Character scalar}. The output format, either
#'   linkmap (\code{"long"}) or matrix (\code{"wide"}). (Default: \code{"long"})
#' 
#' @param mode \code{Character scalar}. The mode of the output, either as
#'   \code{"presence"} or \code{"coverage"} information.
#'   (Default: \code{"presence"})
#' 
#' @param threshold \code{Numeric scalar}. The coverage threshold to infer
#'   presence, between 0 and 1. Only used for complex modules.
#'   (Default: \code{1})
#' 
#' @param verbose \code{Logical scalar}. Should messages be printed in the
#'   console. (Default: \code{TRUE})
#' 
#' @param timeout \code{Numeric scalar}. The timeout for downloading resources.
#'   (Default: \code{1e6})
#' 
#' @param ... Additional arguments.
#' \itemize{
#'     \item \code{batch.size}: \code{Numeric scalar}. The maximum batch size
#'     for SPARQL or API queries. (Default: half the maximum query size)
#'
#'     \item \code{workers}: \code{Numeric scalar}. Number of workers to use,
#'     automatically detected when \code{NULL}. (Default: \code{NULL})
#'     
#'     \item \code{factor}: \code{Numeric scalar}. Number of jobs per worker.
#'     (Default: \code{3})
#' }
#' 
#' @return
#' A two-column data.frame (x-to-y linkmap) or an x-by-y matrix (presence or
#' coverage).
#' 
#' @examples
#' 
#' library(mia)
#' 
#' # Import dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
#' # Load resource graph
#' graph <- ariadne()
#' 
#' # Retrieve taxon names
#' tax.labs <- getTaxonomyLabels(tse, make.unique = FALSE)
#' tax.labs <- sub("^.+:", "", tax.labs)
#' 
#' # Weave path starting from initial values
#' tax2bugsig <- weavePath(graph, taxname ~ bugsig, init = tax.labs)
#' 
#' # Weave 3-rd path
#' tax2bugsig <- weavePath(graph, taxname ~ bugsig, init = tax.labs, k = 3)
#' 
#' Weave path passing through taxid
#' tax2bugsig <- weavePath(graph, taxname ~ taxid ~ bugsig, init = tax.labs)
#' 
NULL


S7::method(weavePath, igraph) <- function(graph, by, k = 1, init = NULL,
    prune = TRUE, output.format = "long", mode = "presence", threshold = 1,
    verbose = TRUE, timeout = 1e6, ...){
    # Check shared character args
    output.format <- match.arg(output.format, c("long", "wide"))
    # Check shared numeric args
    if( !is.numeric(k) || length(k) != 1L || k <= 0 ){
        stop("'k' must be a positive integer.", call. = FALSE)
    }
    if( !is.numeric(timeout) || length(timeout) != 1L || timeout <= 0 ){
        stop("'timeout' must be a positive number", call. = FALSE)
    }
    # Check shared logical args
    if( !is.logical(prune) || length(prune) != 1L ){
        stop("'prune' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !is.logical(verbose) || length(verbose) != 1L ){
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    # Define complex modules
    mod.names <- c("gbm", "gmm")
    # If formula includes complex modules
    if( any(mod.names %in% all.vars(by)) ){
        # Dispatch to method for complex modules
        out <- .weave_complex(
            graph, by, k, init, prune, output.format,
            mode, threshold, verbose, timeout, ...
        )
    }else{
        # Dispatch to regular method
        out <- .weave_path(
            graph, by, k, init, prune, output.format, verbose, timeout, ...
        )
    }
    return(out)
}


#' @importFrom igraph as_data_frame
#' @importFrom MultiFactor MultiFactor weave
.weave_path <- function(graph, by, k, init, prune, output.format, verbose, timeout, ...){
    # Set timeout for downloads
    options(timeout = timeout)
    # Select unique input
    init <- unique(init)
    # Retrieve edges and nodes data
    graph_df <- as_data_frame(graph, what = "both")
    # Draw kth path from source to target
    path_df <- .draw_path(graph, by, k)
    # Initialise list of linkmaps
    linkmaps <- list()
    # Perform step of path
    for( i in seq_len(nrow(path_df)) ){
        # Retrieve step
        g <- path_df[i, ]
        # If pruning is active
        if( g$step != 1L && prune ){
            # Retrieve logs
            logs <- path_df[1:i - 1, ]
            # Update init
            init <- .update_init(g$from, logs, linkmaps)
        }
        # Print step
        if( verbose ) message(g$from, " -(", g$source, ")-> ", g$to)
        # Fetch linkmap
        linkmap <- .fetch_resource(
            graph_df, g$from, g$to, g$source, init, timeout, ...
        )
        # Add to linkmaps
        linkmaps[[paste0(g$from, "2", g$to, ":", g$source)]] <- linkmap
    }
    # Merge linkmaps by colnames
    linkmaps <- .merge_linkmaps(linkmaps)
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    # Extract formula (remove when MF supports multiple ~)
    inner.by <- .extract_formula(by)
    # Weave desired linkmap from MultiFactor
    out <- weave(mf, inner.by)
    # Convert to wide format
    if( output.format == "wide" ){
        # Convert to adjacency matrix
        out <- table(out) == 1
        # Remove dimnames
        names(dimnames(out)) <- NULL
    }
    return(out)
}


.extract_formula <- function(by){
    # Extract by vars
    by.vars <- .formula2list(by)
    from <- paste0(by.vars[[1]], collapse = "+")
    to <- paste0(by.vars[[length(by.vars)]], collapse = "+")
    by <- as.formula(paste0(from, "~", to))
    return(by)
}


.update_init <- function(from, logs, linkmaps){
    
    logs <- logs[logs$to == from, ]
    
    log.names <- paste0(logs$from, "2", logs$to, ":", logs$source)
    logs <- vapply(
        linkmaps[log.names], `[`, from,
        FUN.VALUE = data.frame(1L), USE.NAMES = FALSE
    )
    
    init <- unique(do.call(c, logs))
    return(init)
}


.merge_linkmaps <- function(linkmaps){
    # Create group keys by sorted column names joined by "2"
    groups <- vapply(
        linkmaps,
        function(x) paste(sort(names(x)), collapse = "2"),
        character(1L)
    )
    # Split list by group keys
    linkmaps <- split(linkmaps, groups)
    # Combine elements group-wise
    linkmaps <- lapply(linkmaps, function(x) do.call(rbind, x))
    return(linkmaps)
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
    
    if( g$source %in% c("BugSigDB", "ChocoPhlAn", "GO", "TIGRfams", "WoL")){
        
        cached <- .cache_resource(g$url, g$source, g$from, g$to)
        
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
        
        g$from <- from
        g$to <- to
        
        spec.from <- node_df[node_df$name == g$from, g$source]
        spec.to <- node_df[node_df$name == g$to, g$source]
        
        df <- .queryOTT(spec.from, spec.to, init, timeout, ...)
    
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
