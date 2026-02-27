
#' @name ariadne
#' @rdname ariadne

#' @export
#' @importFrom httr2 request req_perform resp_body_json
#' @importFrom igraph read_graph as_data_frame graph_from_data_frame
#' @importFrom tools R_user_dir
ariadne <- function(versions = NULL){
    
    # Initialise database
    db <- R_user_dir("ariadne", "data")
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
        # Download file
        file_path <- file.path(db, key)
        request(url) |>
            req_perform(path = file_path)
        # Import graph from file
        file_path <- file.path(db, key)
        graph <- read_graph(file_path, format = "gml")
        graph_df <- as_data_frame(graph, what = "both")
        # Store edge and node data
        edge_dfs[[key]] <- graph_df$edges
        node_dfs[[key]] <- graph_df$vertices
    }
    # Build and clean edge data
    edge_df <- do.call(rbind, edge_dfs)
    edge_df$source <- gsub("\\.gml.*$", "", rownames(edge_df))
    rownames(edge_df) <- NULL
    # Build and clean node data
    node_df <- do.call(rbind, node_dfs)
    node_df$source <- gsub("\\.gml.*$", "", rownames(node_df))
    rownames(node_df) <- NULL
    node_df$id <- NULL
    # Widen database-specific names
    node_df <- reshape(
        node_df, idvar = "name", timevar = "source", direction = "wide"
    )
    names(node_df) <- gsub("specific.", "", names(node_df), fixed = TRUE)
    # Retrieve database-specific names
    spec.from <- .generic2specific(edge_df, node_df, "from")
    spec.to <- .generic2specific(edge_df, node_df, "to")
    # Build resource paths
    paths <- mapply(
        FUN = .build_path,
        from = spec.from,
        to = spec.to,
        repo = edge_df$source,
        USE.NAMES = FALSE
    )
    # Add paths to edge data
    edge_df$path <- as.vector(paths)
    # Build final resource graph
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    return(graph)
}


.generic2specific <- function(edges, nodes, what = c("from", "to")){
    # Match
    idx <- match(edges[ , what], nodes$name)
    idy <- match(edges$source, names(nodes))
    
    specific <- nodes[cbind(idx, idy)]
    return(specific)
}

###

# This is temporary while ariadne.db is not online
meta <- list(
    ChocoPhlAn = "https://zenodo.org/records/17100034/files/",
    GM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/",
    GO = "https://current.geneontology.org/ontology/external2go/",
    KEGG = "https://www.genome.jp/kegg/",
    TIGRFAMs = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/",
    UniProt = "https://sparql.uniprot.org/",
    WoL = "https://ftp.microbio.me/pub/wol-20April2021/"
)

meta <- data.frame(
  name = names(meta),
  repo = unlist(meta, use.names = FALSE)
)

###

# Build path to resource
.build_path <- function(from, to, repo){
    
    FUN <- switch(
        repo,
        ChocoPhlAn = function(from, to, repo){
            paste0(repo, "map_", from, "_", to, ".txt.gz")
        },
        GM = function(from, to, repo){
            if( from == "gmm" ){
                file_name <- "GMMs.v1.07.txt"
            }else if( from == "gbm" ){
                file_name <- "GBMs.v1.0.txt"
            }
            paste0(repo, file_name)
        },
        GO = function(from, to, repo){
            paste0(repo, from, "2", to)
        },
        TIGRfams = function(from, to, repo){
          paste0(from, "_", to, "_LINK")
        },
        WoL = function(from, to, repo){
            # Account for exceptions
            prefix <- ifelse(to == "all", "go", to)
            prefix <- ifelse(to == "ko", "kegg", prefix)
            prefix <- ifelse(to == "protein", "metacyc", prefix)
            # Create path
            paste0(repo, "function/", prefix, "/", to, ".map.xz")
        },
        function(from, to, repo) NA
    )
    
    base_url <- meta$repo[match(repo, meta$name)]
    url <- FUN(from, to, base_url)
    return(url)
}
