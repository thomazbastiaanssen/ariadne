#' Build ariadne resource graph
#'
#' @name ariadne
#' 
#' @description
#' \code{ariadne} imports the resource graph hosted in the companion package
#' ariadne.db.
#' 
#' @param versions \code{Character list}. A named list of resource versions to
#' load. Latest versions are used if not specified. (Default: \code{NULL})
#' 
#' @returns An igraph object with the ariadne resource graph.
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link{listResourceVersions}}
#'   \item ariadne.db: \url{https://github.com/Minotau-R/ariadne.db}
#'   \item Zenodo: \url{https://zenodo.org/records/18788725}
#' }
#' 
#' @examples
#' # Import default resource graph
#' graph <- ariadne()
#' 
#' # Specify custom resource versions
#' versions <- list(BugSigDB = "v1.2.2", WoL = "v20April2021")
#' graph <- ariadne(versions = versions)
NULL


#' @export
#' @rdname ariadne
#' @importFrom igraph read_graph as_data_frame graph_from_data_frame
#' @importFrom stats setNames reshape
#' @importFrom BiocParallel bpmapply
#' @importFrom data.table rbindlist dcast
ariadne <- function(versions = NULL){
    # Import version metadata
    meta <- versionMetadata
    # Set default versions
    default_vers <- meta[meta$default, ]
    default_vers <- setNames(default_vers$version, default_vers$source)
    # If custom versions are not specified
    if( is.null(versions) ){
        # Set to default versions
        versions <- default_vers
    }else{
        # Find unspecified resource versions
        to_add <- setdiff(names(default_vers), names(versions))
        # Add default versions for unspecified resources
        versions <- c(versions, default_vers[to_add])
    }
    # Match requested versions to all available versions
    requested <- paste(names(versions), versions)
    available <- paste(meta$source, meta$version)
    idx <- match(requested, available)
    # Check that versions are valid
    if( any(is.na(idx)) ){
        stop(
            "Some 'versions' were not found. If it is a very recent release, ",
            "it may not be registered in the database yet.", call. = FALSE
        )
    }
    # Select desired resource versions
    meta <- meta[idx, ]
    # Build urls to download resource graphs
    urls <- paste0(
        "https://zenodo.org/records/", meta$graph, "/files/", meta$source, ".gml"
    )
    # Fetch individual resource graphs
    graph_dfs <- bpmapply(
        .fetch_graph, meta$key, urls, SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    # Name each graph by corresponding resource
    names(graph_dfs) <- meta$source
    # Build edge data
    edge_df <- graph_dfs |>
        lapply(`[[`, "edges") |>
        rbindlist(idcol = "source", fill = TRUE)
    # Reorder edge columns
    edge_df <- edge_df[ , c("from", "to", "source", "url")]
    # Build node data
    node_df <- graph_dfs |>
        lapply(`[[`, "vertices") |>
        rbindlist(idcol = "source", fill = TRUE)
    # Reduce missing characters to standard NA
    node_df$url[node_df$url == "NA"] <- NA
    # Remove rownames and store node urls
    node_urls <- unique(node_df[!is.na(node_df$url) , c("name", "url")])
    # Widen database-specific names
    node_df <- dcast(node_df, name ~ source, value.var = "specific")
    # Add back node urls
    node_df  <- merge(node_df, node_urls, by = "name", all.x = TRUE)
    # Build final resource graph
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    # Add versions as attribute to graph
    attr(graph, "versions") <- unlist(versions)
    return(graph)
}


.fetch_graph <- function(key, url){
    # Fetch graph from ariadne.db
    graph <- read_graph(url, format = "gml")
    graph_df <- as_data_frame(graph, what = "both")
    # Add version to edge and node urls
    graph_df <- .insert_version(graph_df, key)
    return(graph_df)
}


.insert_version <- function(graph_df, key){
    graph_df <- lapply(graph_df, function(x){
        if(!is.null(x$url)) x$url <- sub("{version}", key, x$url, fixed = TRUE)
        return(x)
    })
    return(graph_df)
}


.generic2specific <- function(edges, nodes, what = c("from", "to")){
    # Match
    idx <- match(edges[[what]], nodes$name)
    idy <- match(edges$source, names(nodes))
    
    specific <- nodes[cbind(idx, idy)]
    return(specific)
}
