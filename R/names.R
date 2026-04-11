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
#' linkNames(graph, "ko", ids = c("K00001", "K00844", "K03455"))
#' 
#' # Search first path including uniref90
#' searchPath(graph, taxname ~ ko, include = "uniref90")
#' 
#' # Search first 5 paths excluding uniref50 and uniref100
#' searchPath(graph, taxname ~ ko, k = 5, exclude = c("uniref50", "uniref100"))
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
    # Get url for name linkmap (or NA for missing)
    url <- node_df$url[node_df$name == x]
    # Map ids and names to one another
    name_links <- .fetch_node(x, url, ids, names)
    # Select ids and names that matched an entry
    name_links <- name_links |>
        .select_matched(x, 1L, ids, verbose) |>
        .select_matched(x, 2L, names, verbose)
    
    # think if there is a better way, maybe earlier step
    linkmap <- bind_rows(
        data.frame(x = if( is.null(ids) ) integer(0L) else ids),
        data.frame(y = if( is.null(names) ) integer(0L) else names)
    )
    
    merge(linkmap)
    # Reset row indices after subsetting
    rownames(name_links) <- NULL
    return(name_links)
})


#' @importFrom data.table fread
#' @importFrom KEGGREST listDatabases keggList keggFind
#' @importFrom stringr str_split
#' @importFrom readr read_lines
#' @importFrom MultiFactor LinkMap
.fetch_node <- function(x, url, ids, names){
    
    # Check nodes where init is necessary
    if( x == "genes" && is.null(ids) && is.null(names) ){
        stop("'ids', 'names' or both must be provided for ", x, ".",
            call. = FALSE)
    }
    
    if( !is.na(url) ){
        
        name_links <- fread(url, header = FALSE)
    
    }else if( x == "bugsig" ){
        
        url <- "https://zenodo.org/records/15272273/files/bugsigdb_signatures_mixed_ncbi.gmt"
        
        name_links <- url |>
            read_lines(skip = 1L) |>
            str_split(fixed("\t")) |>
            vapply(`[`, 1L, FUN.VALUE = character(1L)) |>
            str_split("_", n = 2L, simplify = TRUE) |>
            as.data.frame()
        
        name_links[[1L]] <- sub("bsdb:", "", name_links[[1L]], fixed = TRUE)
    
    }else if( x %in% c(listDatabases(), "ec", "network") ){
        # Use gene ids as input if target is genes
        keys <- if( is.null(ids) ) x else ids
        # Get vector of feature names
        name_vec <- keggList(keys)
        # Maintain only first name (and last for ko)
        to_remove <- ifelse(x == "ko", ".*;", ";.*")
        name_vec <- sub(to_remove, "", name_vec)
        # Convert to linkmap
        name_links <- data.frame(
            x = names(name_vec), y = name_vec, row.names = NULL
        )
    }
    
    name_links <- LinkMap(name_links)
    
    colnames(name_links) <- c(x, paste0(x, ".name"))
    return(name_links)
}


.select_matched <- function(linkmap, x, what, init, verbose){
    
    if( is.null(init) ){
        return(linkmap)
    }
    # Find matches
    idx <- match(init, linkmap[[what]])
    unmatched <- is.na(idx)
    
    if( verbose && any(unmatched) ){
        
        order_fun <- switch(what, identity, rev)
        keys <- order_fun(c("ids", "names"))
        
        warning(sum(unmatched), " ", keys[1L], " for ", x, " ", keys[2L],
            " not found.", call. = FALSE)
    }
    # 
    linkmap <- linkmap[idx[!unmatched], , drop = FALSE]
    return(linkmap)
}
