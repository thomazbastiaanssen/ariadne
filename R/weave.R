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
    if( use.names ) out <- .id2name(graph_df, out, verbose)
    return(out)
})


#' @importFrom data.table fread
#' @importFrom KEGGREST listDatabases keggList
.id2name <- function(graph_df, linkmap, verbose){
    
    node_df <- graph_df$vertices
    target <- colnames(linkmap)[2L]
    
    url <- node_df$url[node_df$name == target]
    
    if( !is.na(url) ){
        
        name.linkmap <- fread(url, header = FALSE)
        
    }else if( target %in% c(listDatabases(), "network") ){
        # Use gene ids as input if target is genes
        ids <- if( target == "genes" ) levels(linkmap[[2L]]) else target
        # Get vector of feature names
        name.vec <- keggList(ids)
        # Maintain only first name
        name.vec <- sub(";.*", "", name.vec)
        # Use gene ids as names if target is genes
        if( target == "genes" ) names(name.vec) <- levels(linkmap[[2L]])
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
    linkmap[paste0(target, ".name")] <- as.factor(name.linkmap[idx, ][[2L]])
    return(linkmap)
}


#' @importFrom KEGGREST keggConv keggLink
#' @importFrom arrow read_parquet open_dataset
#' @importFrom dplyr filter collect
#' @importFrom stats na.omit
#' @importFrom rlang sym
.fetch_resource <- function(graph_df, from, to, repo, init, timeout, ...){
    
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
    
    is_init <- !is.null(init)
    is_source <- edge_df$source == repo
    
    idx <- which(
        edge_df$from == from & edge_df$to == to & is_source
    )
    
    if( length(idx) == 0L ){
        
        idx <- which(
            edge_df$from == to & edge_df$to == from & is_source
        )
    }
    # Only one edge match allowed
    if( length(idx) > 1L ) stop("Multiple matches were found.", call. = FALSE)
    # Retrieve matched edge (only one)
    g <- edge_df[idx, , drop = FALSE]
    # Check edges where init is necessary
    if( !is_init &&
        (g$source == "OTT" || (g$source == "KEGG" && from == "genes")) ){
        stop("'init' must be provided for ", g$from, " queries to ", g$source,
            ".", call. = FALSE)
    }
    # Retrieve specific names for KEGG, OTT and SPARQL queries
    spec.from <- node_df[node_df$name == from, g$source]
    spec.to <- node_df[node_df$name == to, g$source]
    # Retrieve linkmap from corresponding resource
    if( g$source == "KEGG" ){
        # List external databases
        ext <- c("chebi", "geneid", "proteinid", "pubchem", "uniprotkb")
        # Select function based on id types
        kegg_fun <- ifelse(any(c(from, to) %in% ext), keggConv, keggLink)
        # Use initial values as input for genes db
        orig <- if( "genes" %in% c(from, to) ) init else spec.from
        # Add prefix to external from ids
        if( is_init && from %in% ext ) orig <- paste0(spec.from, ":", orig)
        # Send query to keggLink
        kegg_link <- kegg_fun(spec.to, orig)
        # Convert to data.frame
        df <- data.frame(x = names(kegg_link), y = kegg_link, row.names = NULL)
        # Strip db prefix except for genes db
        if( from != "genes") df$x <- sub("^[^:]*:", "", df$x)
        if( to != "genes" ) df$y <- sub("^[^:]*:", "", df$y)
        # Use initial values to filter output
        if( is_init ) df <- df[df[[1L]] %in% init, , drop = FALSE]
    # Query Open Tree Taxonomy API
    }else if( g$source == "OTT" ){
        # Send OTT query (init must be vector)
        df <- .queryOTT(spec.from, spec.to, init, timeout, ...)
    # Query SPARQL endpoint
    }else if( g$source %in% c("Rhea", "UniProt") ){
        # Add special IRI prefixes
        if( is_init ) init <- .add_iri(init, g$source, from)
        # Query SPARQL endpoint
        df <- .querySPARQL(spec.from, spec.to, g$source, init, timeout, ...)
        # Filter special cases
        print(head(df))
        if( spec.to %in% c("uniref", "BioCyc") ){
            df <- df[grepl(to, df[[2L]], ignore.case = TRUE), ]
            rownames(df) <- NULL
        }
        print(head(df))
        # Strip special IRI prefixes
        df[[1L]] <- .strip_iri(df[[1L]], from)
        df[[2L]] <- .strip_iri(df[[2L]], to)
    # Fetch linkmap from file
    }else{
        # Get file path to cached resource
        cached <- .cache_resource(g$url, g$source, g$from, g$to)
        df <- read_parquet(cached)
        # If initial values are given
        if( is_init ){
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
        # Replace colnames with specifics
        from <- g$from
        to <- g$to
    }
    # Check that result is not empty
    if( nrow(df) == 0L ) stop("Bindings depleted.", call. = FALSE)
    # Add edge names
    colnames(df) <- c(from, to)
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
