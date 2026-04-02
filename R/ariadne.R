#' Build ariadne resource graph
#'
#' @name ariadne
#' @rdname ariadne
#' 
#' @description
#' \code{ariadne} imports the resource graph hosted in ariadne.db companion data
#' package.
#' 
#' @param versions \code{Character list}. A named list of resource versions to
#' load. Latest versions are used if not specified. (Default: \code{NULL})
#' 
#' @returns
#' An igraph object with the ariadne resource graph.
#' 
#' @examples
#' 
#' # Import ariadne resource graph
#' graph <- ariadne()
#' 
#' # Specify resource versions
#' graph <- ariadne(versions = list(BugSigDB = "v1.2.0", WoL = "v2.0"))
#' 
#' @seealso
#' 
#' ariadne.db: \url{https://github.com/Minotau-R/ariadne.db}
#' Zenodo: \url{https://zenodo.org/records/18788725}
#' 
NULL

#' @export
#' @importFrom httr2 request req_perform resp_body_json
#' @importFrom igraph read_graph as_data_frame graph_from_data_frame
#' @importFrom dplyr bind_rows
#' @importFrom stats reshape
ariadne <- function(versions = NULL){
    # Define database url
    url <- "https://zenodo.org/api/records/18788725"
    # Send request to database
    resp <- url |>
        request() |>
        req_perform()
    # Check if request was successful
    if( resp$status_code != 200 ){
        stop("Failed to retrieve ariadne.db.", call. = FALSE)
    }
    # Parse JSON content
    record_json <- resp_body_json(resp)
    # Extract file download URLs and filenames
    files <- record_json$files
    keys <- vapply(files, `[[`, "key", FUN.VALUE = character(1L))
    urls <- vapply(files, function(x) x$links$self, character(1L))
    # Initialise edge and node data
    edge_dfs <- list()
    node_dfs <- list()
    # For each resource graph
    for( i in seq_along(keys) ){
        # Retrieve current key and url
        key <- keys[i]
        url <- urls[i]
        # Fetch graph from ariadne.db
        graph <- read_graph(url, format = "gml")
        graph_df <- as_data_frame(graph, what = "both")
        # Store edge and node data
        edge_dfs[[key]] <- graph_df$edges
        node_dfs[[key]] <- graph_df$vertices
    }
    # Build and clean edge data
    edge_df <- bind_rows(edge_dfs, .id = "source")
    edge_df$source <- sub(".gml", "", edge_df$source, fixed = TRUE)
    edge_df <- edge_df[ , c("from", "to", "source", "url")]
    # Build and clean node data
    node_df <- bind_rows(node_dfs, .id = "source")
    node_df$source <- sub(".gml", "", node_df$source, fixed = TRUE)
    node_df$url[node_df$url == "NA"] <- NA
    # Remove rownames and store node urls
    rownames(node_df) <- NULL
    node_urls <- unique(node_df[c("name", "url")])
    # Widen database-specific names
    node_df <- reshape(
        node_df, idvar = "name", timevar = "source",
        direction = "wide", drop = c("id", "url")
    )
    # Clean colnames and add back node urls
    names(node_df) <- sub("specific.", "", names(node_df), fixed = TRUE)
    node_df  <- merge(node_df, node_urls, by = "name", all.x = TRUE)
    # Build final resource graph
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    return(graph)
}


.generic2specific <- function(edges, nodes, what = c("from", "to")){
    # Match
    idx <- match(edges[[what]], nodes$name)
    idy <- match(edges$source, names(nodes))
    
    specific <- nodes[cbind(idx, idy)]
    return(specific)
}
