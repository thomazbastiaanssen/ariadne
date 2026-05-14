#' Weave paths between resources
#' 
#' @name weavePath
#' @aliases weaveComplex
#' 
#' @description
#' \code{weavePath} and \code{weaveComplex} build the selected path through the
#' resource graph from the origin to the target features, fetching and combining
#' the necessary resources from several databases.
#' 
#' While \code{weavePath} returns a linkmap with all simple links,
#' \code{weaveComplex} yields links only for features mapped above a certain
#' coverage threshold, which is useful for pathways or functional modules made
#' of several indispensable components.
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
#' @param res.name \code{Character vector}. Names of resources to include in
#'   the graph. (Default: \code{NULL})
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
#' @param threshold \code{Numeric scalar}. Only for \code{weaveComplex}. The
#'   coverage threshold above which to include links, between 0 and 1. If
#'   \code{NULL}, all non-zero links are returned. (Default: \code{NULL})
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
#' @returns
#' A two-column data.frame (x-to-y linkmap), with optional names as a third
#' column.
#' 
#' @examples
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
#' dis2gmm <- weavePath(graph, kegg_disease ~ gmm)
#' 
#' # Weave complex path from KEGG diseases to gut metabolic modules
#' dis2gmm <- weaveComplex(graph, kegg_disease ~ gmm)
#' 
#' # Specify coverage threshold
#' dis2gmm <- weaveComplex(graph, kegg_disease ~ gmm, threshold = 0.8)
NULL


#' @export
#' @rdname weavePath
#' @importFrom stats as.formula
#' @importFrom MultiFactor MultiFactor weave
setMethod("weavePath", signature = c(graph = "igraph"),
    function(graph, by, k = 1, include = NULL, exclude = NULL, res.name = NULL,
    init = NULL, prune = TRUE, use.names = TRUE, verbose = TRUE,
    timeout = 1e6, ...){
    # Build MultiFactor from path linkmaps
    mf <- .build_path_mf(
        graph, by, k, include, exclude, res.name,
        init, prune, prune, verbose, timeout, ...
    )
    # Weave desired linkmap from MultiFactor
    out <- weave(mf, by) |> as.data.frame()
    # Add feature names
    if( use.names ){
        target <- colnames(out)[2L]
        name_links <- linkNames(graph, target, out[[2L]], verbose = verbose)
        out[paste0(target, ".name")] <- as.factor(name_links[[2L]])
    }
    return(out)
})


.build_path_mf <- function(graph, by, k, include, exclude, res.name, init,
    prune, prune.last, verbose, timeout, ...){
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
    # Draw kth path from source to target
    path_df <- .draw_path(graph, path_by, k, include, exclude, res.name)
    # Add edges metadata
    path_df <- .add_edge_metadata(path_df, graph, internal = TRUE)
    # Create pruning instructions
    prune_vec <- c(rep(prune, max(0, nrow(path_df) - 2)), prune.last, FALSE)
    # Perform step of path
    for( i in seq_len(nrow(path_df)) ){
        # Retrieve step
        g <- path_df[i, , drop = FALSE]
        # Print step
        if( verbose ) message(g$initFrom, " -(", g$source, ")-> ", g$initTo)
        # Fetch linkmap
        linkmap <- .fetch_edge(g, init, timeout, ...)
        # Add to linkmaps
        linkmaps[[paste0(g$from, "2", g$to)]] <- linkmap
        # Update init
        init <- if( prune_vec[i] ) levels(linkmap[[g$initTo]]) else NULL
    }
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    return(mf)
}


#' @importFrom KEGGREST keggConv keggLink
#' @importFrom arrow read_parquet open_dataset
#' @importFrom dplyr filter collect
#' @importFrom MultiFactor LinkMap
#' @importFrom rlang sym
.fetch_edge <- function(g, init, timeout, ...){
    # Check if init exists
    is_init <- !is.null(init)
    # Check edges where init is necessary
    if( !is_init &&
        (g$source == "OTT" || (g$source == "KEGG" && g$from == "kegg_genes")) ){
        stop("'init' must be provided for ", g$from, " queries to ", g$source,
            ".", call. = FALSE)
    }
    # Retrieve linkmap from corresponding resource
    if( g$source == "KEGG" ){
        # List external databases
        ext <- c("chebi", "geneid", "proteinid", "pubchem", "uniprotkb")
        # Select function based on id types
        kegg_fun <- ifelse(any(c(g$from, g$to) %in% ext), keggConv, keggLink)
        # Use initial values as input for genes db
        orig <- if( "kegg_genes" %in% c(g$from, g$to) ) init else g$specFrom
        # Add prefix to external from ids
        if( is_init && g$from %in% ext ) orig <- paste0(g$specFrom, ":", orig)
        # Send query to keggLink
        kegg_link <- kegg_fun(g$specTo, orig)
        # Convert to data.frame
        df <- data.frame(x = names(kegg_link), y = kegg_link, row.names = NULL)
        # Strip db prefix except for genes db
        if( g$from != "kegg_genes") df$x <- sub("^[^:]*:", "", df$x)
        if( g$to != "kegg_genes" ) df$y <- sub("^[^:]*:", "", df$y)
        # Use initial values to filter output
        if( is_init ) df <- df[df[[1L]] %in% init, , drop = FALSE]
    # Query Open Tree Taxonomy API
    }else if( g$source == "OTT" ){
        # Send OTT query (init must be vector)
        df <- .queryOTT(g$specFrom, g$specTo, init, timeout, ...)
    # Query SPARQL endpoint
    }else if( g$source %in% c("Rhea", "UniProt") ){
        # Add special IRI prefixes
        if( is_init ) init <- .add_iri(init, g$source, g$from)
        # Select similarity level for uniref clusters
        uniref.identity <- switch(g$to, uniref50 = 0.5, uniref90 = 0.9, NULL)
        # Query SPARQL endpoint
        df <- .querySPARQL(
            g$specFrom, g$specTo, g$source, init, uniref.identity, timeout, ...
        )
        # Filter special cases
        if( g$specTo == "BioCyc" ){
            # Get key (metacyc or ecocyc)
            biocyc_key <- sub("_.+$", "", g$to)
            # Filter by key
            df <- df[grepl(biocyc_key, df[[2L]], ignore.case = TRUE), ]
            rownames(df) <- NULL
        }
        # Strip special IRI prefixes
        df[[1L]] <- .strip_iri(df[[1L]], g$specFrom)
        df[[2L]] <- .strip_iri(df[[2L]], g$specTo)
    # Fetch linkmap from file
    }else{
        # Get file path to cached resource
        cached <- .cache_resource(g$url, g$source, g$from, g$to)
        # If initial values are given
        if( is_init ){
            # Filter linkmap before importing
            df <- cached |>
                open_dataset() |>
                filter(!!sym(g$initFrom) %in% init) |>
                collect() |>
                as.data.frame()
        }else{
            # Read linkmap from parquet
            df <- read_parquet(cached)
        }
        # Swap columns if direction does not equal file order 
        if( g$from != g$initFrom ) df <- rev(df)
    }
    # Check that result is not empty
    if( nrow(df) == 0L ) stop("Bindings depleted.", call. = FALSE)
    # Add edge names
    colnames(df) <- c(g$initFrom, g$initTo)
    # Convert characters to factors
    df <- LinkMap(df)
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
    
    if( !grepl("^genes$|^uniref", name) ){
        y <- sub(paste0("^", name, "[:_]"), "", y, ignore.case = TRUE)
    }
    
    return(y)
}
