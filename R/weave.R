#' Weave paths between resources
#' 
#' @description
#' \code{weavePath} and \code{weaveComplex} bridge the path between resources,
#' fetching and combining the necessary data from the ariadne database.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to weave.
#' 
#' @param k \code{Numeric scalar}. The kth shortest path to weave.
#'   (Default: \code{1})
#' 
#' @param include \code{Character vector}. Nodes to cross in the path.
#'   (Default: \code{NULL})
#' 
#' @param exclude \code{Character vector}. Nodes to avoid in the path.
#'   (Default: \code{NULL})
#' 
#' @param init \code{Character vector}. Initial values to prune the first
#'   mapping step. (Default: \code{NULL})
#' 
#' @param prune \code{Logical scalar}. Should the values from each step be used
#'   prune the following step. (Default: \code{TRUE})
#' 
#' @param use.names \code{Logical scalar}. Should feature names be used in the
#'   output instead of feature identifiers. either (Default: \code{TRUE})
#' 
#' @param mode \code{Character scalar}. The mode of the output, either as
#'   \code{"presence"} or \code{"coverage"} information. Only for
#'   \code{weaveComplex}. (Default: \code{"presence"})
#' 
#' @param threshold \code{Numeric scalar}. The coverage threshold to infer
#'   presence, between 0 and 1. Only for \code{weaveComplex}.
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
#'   \item \code{batch.size}: \code{Numeric scalar}. The maximum batch size
#'   for SPARQL or API queries. (Default: half the maximum query size)
#'
#'   \item \code{workers}: \code{Numeric scalar}. Number of workers to use,
#'   automatically detected when \code{NULL}. (Default: \code{NULL})
#'     
#'   \item \code{factor}: \code{Numeric scalar}. Number of jobs per worker.
#'   (Default: \code{3})
#' }
#' 
#' @returns A two-column data.frame (x-to-y linkmap).
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
#' # Weave path including taxid
#' tax2bugsig <- weavePath(
#'     graph, taxname ~ bugsig, include = "taxid", init = tax.labs
#' )
#' 
#' # Weave simple path from KEGG diseases to gut metabolic modules
#' dis2gmm <- weavePath(graph, disease ~ gmm)
#' 
#' # Weave complex path from KEGG diseases to gut metabolic modules
#' dis2gmm <- weaveComplex(graph, disease ~ gmm, threshold = 0.8)
#' 
#' # Obtain results in terms of coverage
#' dis2gmm <- weaveComplex(graph, disease ~ gmm, mode = "coverage")
#' 
#' @name weavePath
#' @aliases weaveComplex
NULL


#' @export
#' @rdname weavePath
#' @importFrom igraph as_data_frame
#' @importFrom stats as.formula
#' @importFrom MultiFactor MultiFactor weave
setMethod("weavePath", signature = c(graph = "igraph"),
    function(graph, by, k = 1, include = NULL, exclude = NULL, init = NULL,
    prune = TRUE, use.names = TRUE, verbose = TRUE, timeout = 1e6, ...){
    
    # Check shared numeric args
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
    # Set timeout for downloads
    options(timeout = timeout)
    # Initialise list of linkmaps
    linkmaps <- list()
    # For stratified input
    if( is.data.frame(init) && ncol(init) == 2L ){
        # Retrieve init variable names
        init_vars <- colnames(init)
        # Remove first step in the path
        path_by <- as.formula(paste0(init_vars[2L], "~", all.vars(by)[2L]))
        # Print stratification
        if( verbose ) message(init_vars[1L], " stratified by ", init_vars[2L])
        # Add init linkmap to linkmaps
        linkmaps[["init"]] <- init
        # Extract initial values for second step
        init <- init[[2L]]
    }else{
        path_by <- by
    }
    # Select unique input
    init <- unique(init)
    # Retrieve edges and nodes data
    graph_df <- as_data_frame(graph, what = "both")
    # Draw kth path from source to target
    path_df <- .draw_path(graph, path_by, k, include, exclude)
    # Perform step of path
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
    out <- weave(mf, by)
    # Add feature names
    out <- if( use.names ) .id2name(graph_df, out, verbose) else out
    return(out)
})


#' @importFrom utils read.table
#' @importFrom KEGGREST keggList
.id2name <- function(graph_df, linkmap, verbose){
    
    node_df <- graph_df$vertices
    target <- colnames(linkmap)[2L]
    
    url <- node_df$url[node_df$name == target]
    
    if( !is.na(url) ){
        
        name.linkmap <- read.table(url, sep = "\t")
        
    }else if( target %in% node_df$KEGG[!is.na(node_df$KEGG)] ){
        # Get vector of feature names
        name.vec <- keggList(target)
        # Maintain only first name
        name.vec <- sub(";.*", "", name.vec)
        # Convert to linkmap
        name.linkmap <- data.frame(
            x = names(name.vec), y = name.vec, row.names = NULL
        )
    }else{
        return(linkmap)
    }
    # Find matches
    idx <- match(linkmap[[2L]], name.linkmap[[1L]])
    
    unmatched <- is.na(idx)
    
    if( verbose && any(unmatched) ){
        warning("Names for ", sum(unmatched), " ", target, " ids not found.",
            call. = FALSE)
    }
    # Map ids to names
    linkmap[paste0(target, ".name")] <- as.factor(name.linkmap[idx, 2L])
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
    
    is.init <- !is.null(init)
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
    
    if( g$source %in% c("BugSigDB", "ChocoPhlAn", "GM", "GO", "TIGRFAMs", "WoL")){
        # Get file path to cached resource
        cached <- .cache_resource(g$url, g$source, g$from, g$to)
        # If initial values are given
        if( is.init ){
            # Filter linkmap before importing
            df <- cached |>
                open_dataset() |>
                filter(!!sym(from) %in% init) |>
                collect() |>
                as.data.frame()
        }else{
            # Read linkmap from parquet
            df <- read_parquet(cached)
        }
        
    }else if( g$source == "KEGG" ){
    
        kegg.link <- keggLink(g$from, g$to)
        
        df <- data.frame(
            x = gsub("^[^:]*:", "", kegg.link),
            y = gsub("^[^:]*:", "", names(kegg.link))
        )
        
        if( is.init ){
            df <- df[df[[1L]] %in% init, , drop = FALSE]
        }
        
    }else if( g$source == "OTT" ){
        
        g$from <- from
        g$to <- to
        
        spec.from <- node_df[node_df$name == g$from, g$source]
        spec.to <- node_df[node_df$name == g$to, g$source]
        
        df <- .queryOTT(spec.from, spec.to, init, timeout, ...)
    
    }else if( g$source %in% c("Rhea", "UniProt") ){
        # Variable original order matters
        g$from <- from
        g$to <- to
        # Add special IRI prefixes
        init <- if( is.init ) .add_iri(init, g$source, g$from) else init
        # Use specific names for SPARQL queries
        spec.from <- node_df[node_df$name == g$from, g$source]
        spec.to <- node_df[node_df$name == g$to, g$source]
        # Query SPARQL endpoint
        df <- .querySPARQL(spec.from, spec.to, g$source, init, timeout, ...)
        # Filter special cases
        if( spec.to %in% c("uniref", "BioCyc") ){
            df <- df[grepl(g$to, df[[spec.to]], ignore.case = TRUE), ]
            rownames(df) <- NULL
        }
        # Strip special IRI prefixes
        df[[spec.from]] <- .strip_iri(df[[spec.from]], g$from)
        df[[spec.to]] <- .strip_iri(df[[spec.to]], g$to)
    }
    # Check that result is not empty
    if( nrow(df) == 0L ){
        stop("Bindings depleted.", call. = FALSE)
    }
    # Add edge names
    colnames(df) <- c(g$from, g$to)
    return(df)
}


.add_iri <- function(x, repo, from){
    if( from %in% c("ecocyc", "metacyc") ){
        prefix <- ifelse(
            repo == "UniProt",
            switch(from, metacyc = "MetaCyc", ecocyc = "EcoCyc"),
            toupper(from)
        )
        x <- paste0(prefix, ":", x)
    }else if( from %in% c("chebi", "go") ){
        x <- paste0(toupper(from), "_", x)
    }
    return(x)
}

.strip_iri <- function(x, name){
    
    y <- sub("http.+/", "", x)
    
    if( !grepl("uniref", name, fixed = TRUE) ){
        y <- sub(paste0("^", name, "[:_]"), "", y, ignore.case = TRUE)
    }
    
    return(y)
}
