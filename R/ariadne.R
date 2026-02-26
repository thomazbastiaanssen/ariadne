
#' @name ariadne
#' @rdname ariadne

#' @importFrom igraph graph_from_data_frame
# Fix method def
S7::method(ariadne) <- function(versions = NULL){

    record_id <- "18788726"
    url <- paste0("https://zenodo.org/api/records/", record_id)
    
    resp <- url |>
        request() |>
        req_perform()
    
    # Check if request was successful
    if( resp$status_code != 200 ){
        stop("Failed to retrieve ariadne.db.", call. = FALSE)
    }
    
    # Parse JSON content
    record_json <- resp |> resp_body_json()
    
    # Extract file download URLs and filenames
    files <- record_json$files
    keys <- vapply(files, function(x) gsub(".zip", "", x$key), character(1L))
    
    # Folder to save files
    db <- R_user_dir("ariadne", "data")
    
    if (!dir.exists(db)) dir.create(db, recursive = TRUE)
    
    message("Fetching database: ", paste0(keys, collapse = ", "))
    
    edge_dfs <- list()
    node_dfs <- list()
    
    # Download each file
    for (key in keys) {
      
        file_name <- paste0(key, ".zip")
        file_path <- file.path(db, file_name)
        dir_path <- file_path_sans_ext(file_path)
        # Download file
        request(file.path(url, "files", file_name, "content")) |>
            req_perform(path = file_path)
        # Unzip contents
        unzip(file_path, exdir = dir_path)
        # Remove zip file
        file.remove(file_path)
        
        edge_df <- read.table(
            file.path(dir_path, "edges.tsv"),
            sep = "\t", header = TRUE
        )
        node_df <- read.table(
            file.path(dir_path, "nodes.tsv"),
            sep = "\t", header = TRUE
        )
        
        edge_dfs[[key]] <- edge_df
        node_dfs[[key]] <- node_df
    }
    
    edge_df <- do.call(rbind, lapply(edge_dfs, `[`))
    edge_df$source <- gsub("\\.[0-9]*$", "", rownames(edge_df))
    rownames(edge_df) <- NULL
    
    node_df <- do.call(rbind, lapply(node_dfs, `[`))
    node_df$source <- gsub("\\.[0-9]*$", "", rownames(node_df))
    rownames(node_df) <- NULL
    
    node_df <- reshape(
        node_df, idvar = "generic", timevar = "source", direction = "wide"
    )
    
    names(node_df) <- gsub("specific.", "", names(node_df), fixed = TRUE)
    
    spec.from <- .generic2specific(edge_df, node_df, "from")
    spec.to <- .generic2specific(edge_df, node_df, "to")
    
    paths <- mapply(
        FUN = .build_path,
        from = spec.from,
        to = spec.to,
        repo = edge_df$source,
        USE.NAMES = FALSE
    )
    
    edge_df$path <- as.vector(paths)
    
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    return(graph)
}

.generic2specific <- function(edges, nodes, what = c("from", "to")){
    # Match
    idx <- match(edges[ , what], nodes$generic)
    idy <- match(edges$source, names(nodes))
    
    specific <- nodes[cbind(idx, idy)]
    return(specific)
}


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
