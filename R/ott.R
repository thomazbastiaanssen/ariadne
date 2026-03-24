
#' @importFrom BiocParallel bplapply
.queryOTT <- function(from, to, init, timeout, ...){
    # Remove rank prefixes
    x <- gsub("^[a-z]__", "", init)
    # If source is taxname
    if( from == "taxname" ){
        # Search with TNRS API
        y <- .queryTNRS(x, to, ...)
    }
    # If source is ott
    else if( from == "ott" ){
        # Convert to numeric
        x <- as.numeric(x)
        # Search with taxon_info API
        y <- bplapply(
            x, .queryTaxonInfo, to = to, prefix = "ott_id", timeout = timeout
        )
    # If source is external taxon id
    }else if( from %in% c("ncbi", "gbif", "worms", "if", "irmng") ){
        # Add prefix
        x <- paste0(from, ":", x)
        # Search with taxon_info API
        y <- bplapply(
            x, .queryTaxonInfo, to = to, prefix = "source_id", timeout = timeout
        )
    }else{
        stop("'from' is not valid.", call. = FALSE)
    }
    # Fill empty results
    y[lengths(y) == 0L] <- NA
    # Create linkmap
    linkmap <- data.frame(
        x = rep(init, lengths(y, use.names = FALSE)),
        y = unlist(y, recursive = TRUE, use.names = FALSE)
    )
    # Omit rows with missing values
    linkmap <- na.omit(linkmap)
    # Strip source prefix
    linkmap$y <- sub("^.+:", "", linkmap$y)
    return(linkmap)
}


#' @importFrom httr2 request req_method req_body_json req_perform resp_body_json
.queryTaxonInfo <- function(x, to, prefix, timeout){
    
    url <- "https://api.opentreeoflife.org/v3/taxonomy/taxon_info"
    
    body <- list()
    body[[prefix]] <- x
    
    resp <- request(url) |>
        req_method("POST") |>
        req_body_json(body) |>
        req_timeout(timeout) |>
        req_perform()
    
    resp <- resp_body_json(resp)
    
    # Handle missing or empty responses
    if( is.null(resp) || length(resp) == 0L ){
        return(NA_character_)
    }
    
    if( to == "ott" ){
        y <- resp$ott_id
    }else if( to == "taxname" ){
        tax_rank <- substr(resp$rank, 1, 1)
        tax_name <- resp$unique_name
        y <- paste0(tax_rank, "__", tax_name)
    }else{
        tax_sources <- unlist(resp$tax_sources)
        y <- tax_sources[grepl(to, tax_sources, fixed = TRUE)]
    }
    return(y)
}

#' @importFrom BiocParallel bplapply
.queryTNRS <- function(x, to, batch.size = 1000, workers = NULL, factor = 3){
    
    ranges <- .get_batches(x, batch.size, workers, factor)
    
    out.list <- bplapply(
       ranges, function(i) .subqueryTNRS(x[i[1]:i[2]], to = to)
    )
    
    out <- do.call(c, out.list)
    return(out)
}

#' @importFrom rotl tnrs_match_names tax_sources
.subqueryTNRS <- function(x, to){
    # Match input names
    res <- tnrs_match_names(x)
    # Based on target
    if( to == "ott" ){
        # Extract ott ids
        y <- res$ott_id
    }else{
        # Find ids for each input
        y <- lapply(tax_sources(res), function(x) x[grepl(to, x, fixed = TRUE)])
        # Match ids to input names
        y <- y[res$unique_name]
    }
    return(y)
}
