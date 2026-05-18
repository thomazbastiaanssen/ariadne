#' Link ids to names and vice versa
#' 
#' @name linkNames
#' 
#' @description
#' \code{linkNames} retrieves the corresponding names for a set of ids and
#' vice versa. It can be used to link feature names to ids that can be directly
#' passed to `weavePath` and `weaveComplex`.
#' 
#' @param graph An igraph object.
#' 
#' @param x \code{Character scalar}. A feature type present as a node in
#'   \code{graph}.
#' 
#' @param ids \code{Character vector}. A set of ids to link to their
#'   corresponding names. (Default: \code{NULL})
#' 
#' @param names \code{Character vector}. A set of names to link to their
#'   corresponding ids. (Default: \code{NULL})
#' 
#' @param verbose \code{Logical scalar}. Should messages be printed in the
#'   console. (Default: \code{TRUE})
#' 
#' @param ... Unused.
#' 
#' @returns
#' A \code{\link[MultiFactor:LinkMap]{LinkMap}} with ids and names in the first
#' and second column, respectively.
#' 
#' @examples
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Fetch names for all Gut Metabolic modules
#' linkNames(graph, "gmm")
#' 
#' # Fetch names for a set of KO ids
#' linkNames(graph, "ko", ids = c("K00001", "K00844", "K03455"))
#' 
#' # Fetch names for a set of EC names
#' linkNames(graph, "ec", names = c("alcohol dehydrogenase", "alditol oxidase"))
NULL


#' @export
#' @rdname linkNames
#' @importFrom igraph as_data_frame
setMethod("linkNames", signature = c(graph = "igraph"),
    function(graph, x, ids = NULL, names = NULL, verbose = TRUE){
    # Retrieve nodes data
    node_df <- as_data_frame(graph, what = "vertices")
    # Check that node exists
    if( !x %in% node_df$name ){
        stop("'x' must be in 'graph'.", call. = FALSE)
    }
    # Check input keys
    if( !is.null(ids) && !is.null(names) ){
        stop("Either 'ids' or 'names' can be specified.", call. = FALSE)
    }
    # Retrieve corresponding node from graph
    g <- node_df[node_df$name == x, , drop = FALSE]
    # Get specific name of x for source database
    g$spec <- .generic2specific(g, node_df, "name")
    # Map ids and names to one another
    name_links <- .fetch_node(g, ids)
    # Return empty object if no names are available
    if( is.null(name_links) ){
        return(NULL)
    }
    # Select ids and names that matched an entry
    name_links <- name_links |>
        .match_key2val(x, 1L, ids, verbose) |>
        .match_key2val(x, 2L, names, verbose)
    # Use keywords as column names
    colnames(name_links) <- c(x, paste0(x, ".name"))
    # Convert to data.frame
    name_links <- as.data.frame(name_links)
    return(name_links)
})


#' @importFrom KEGGREST listDatabases keggList
#' @importFrom arrow read_parquet open_dataset
#' @importFrom dplyr select filter collect
#' @importFrom stringr str_split fixed
#' @importFrom tidyselect all_of
#' @importFrom readr read_lines
#' @importFrom data.table fread
#' @importFrom rlang sym
.fetch_node <- function(g, ids){
    # Check if init ids exist
    is_init <- !is.null(ids)
    # Check nodes where init is necessary
    if( g$name == "kegg_genes" && !is_init ){
        stop("Only searches with 'ids' are currently supported for ", g$name,
            ".", call. = FALSE)
    }
    
    if( g$name == "bugsig" ){
        
        url <- "https://zenodo.org/records/15272273/files/bugsigdb_signatures_mixed_ncbi.gmt"
        
        name_links <- url |>
            read_lines(skip = 1L) |>
            str_split(fixed("\t")) |>
            vapply(`[`, 1L, FUN.VALUE = character(1L)) |>
            str_split(fixed("_"), n = 2L, simplify = TRUE) |>
            as.data.frame()
        
        name_links[[1L]] <- sub("bsdb:", "", name_links[[1L]], fixed = TRUE)
        name_links[[2L]] <- sub("^.+:", "", name_links[[2L]])
    
    }else if( g$name %in% c("gmm", "gbm") ){
      
        name_links <- fread(g$url, header = FALSE, showProgress = FALSE)
    
    }else if( !is.na(g$url) ){
        # Build name colname
        name_col <-  paste0(g$spec, "_name")
        if( g$name == "msig" ) name_col <- sub("_id", "", name_col, fixed = TRUE)
        # Get file path to cached resource
        cached <- .cache_resource(g$url, g$source, g$spec, name_col)
        # If initial values are given
        if( is_init ){
            # Filter linkmap before importing
            name_links <- cached |>
                open_dataset() |>
                dplyr::select(all_of(c(g$spec, name_col))) |>
                filter(!!sym(g$spec) %in% ids) |>
                collect() |>
                as.data.frame()
        }else{
            # Read linkmap from parquet
            name_links <- read_parquet(
                cached, col_select = all_of(c(g$spec, name_col))
            )
        }
    }else if( g$KEGG %in% c(listDatabases(), "ec", "network") ){
        # Use ids as input if specified
        init <- if( is_init ) unique(ids) else g$KEGG
        # For many ids, global search is faster
        if( length(init) > 50 ) init <- g$KEGG
        # Get vector of feature names
        name_vec <- keggList(init)
        # Keep only first name (and last for ko)
        to_remove <- ifelse(g$KEGG == "ko", ".*;", ";.*")
        name_vec <- sub(to_remove, "", name_vec)
        # Convert to linkmap
        name_links <- data.frame(
            x = names(name_vec), y = name_vec, row.names = NULL
        )
    }else if( g$name == "chebi" ){
        query <- "
            PREFIX up: <http://purl.uniprot.org/core/>
            SELECT DISTINCT ?chebi ?name
            WHERE {?chebi up:name ?name}
        "
        name_links <- .sendSPARQL(query, "Rhea", 1e6)
        name_links[[1L]] <- sub(
            "http://purl.obolibrary.org/obo/CHEBI_", "",
            name_links[[1L]], fixed = TRUE
        )
    }else{
        return(NULL)
    }
    # Add placeholders to column names
    colnames(name_links) <- c("ids", "names")
    return(name_links)
}


.match_key2val <- function(linkmap, x, what, init, verbose){
    # Return original linkmap if no keyword is defined
    if( is.null(init) ){
        return(linkmap)
    }
    # Retrieve keywords ("ids" and "names")
    keys <- colnames(linkmap)
    # Define function to order keywords
    order_fun <- switch(what, identity, rev)
    # Order keywords
    keys <- order_fun(keys)
    # Find matched and total unmatched keys
    idx <- match(init, linkmap[[what]])
    total_unmatched <- sum(is.na(idx))
    # If there are unmatched keys
    if( verbose && total_unmatched != 0L ){
        # Warn about unmatched keys
        warning(keys[2L], " for ", total_unmatched, " ",  x, " ", keys[1L],
            " not found.", call. = FALSE)
    }
    # Create id2name linkmap
    linkmap <- data.frame(
        x = init, y = linkmap[idx, ][[keys[2L]]], row.names = NULL
    )
    # Reorder columns based on input keyword
    linkmap <- order_fun(linkmap)
    return(linkmap)
}
