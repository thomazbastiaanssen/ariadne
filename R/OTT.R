
# ott, ncbi, gbif, worms, if, irmng, taxname
# silva

# Query from ncbi ids
#ncbi_ids <- c(562, 1423, 1280)
#ott_ids <- .queryOTT(ncbi_ids, from = "ncbi", to = "ott")
#tax_names <- .queryOTT(ncbi_ids, from = "ncbi", to = "taxname")

# Qeury from ott ids
#ott_ids <- c(474506, 1084928, 1090496)
#ncbi_ids <- .queryOTT(ott_ids, from = "ott", to = "ncbi")
#tax_names <- .queryOTT(ott_ids, from = "ott", to = "taxname")

# Query silva ids
#silva_ids <- .queryOTT(ott_ids, from = "ott", to = "silva")
# silva_ids <- .queryOTT(ncbi_ids, from = "ncbi", to = "silva")

# Query from taxnames
#names <- c("s__Escherichia coli", "s__Bacillus subtilis", "s__Staphylococcus aureus")
#ott_ids <- .queryOTT(names, from = "taxname", to = "ott")
#silva_ids <- .queryOTT(names, from = "taxname", to = "silva")


#' @importFrom BiocParallel bplapply
.queryOTT <- function(x, from, to, timeout, ...){
    #
    x <- gsub("^.__", "", x)
    # 
    if( from == "taxname" ){
        y <- .queryTNRS(x, to, timeout, ...)
    }
    else if( from == "ott" ){
        x <- as.numeric(x)
        y <- unlist(bplapply(
            x, .queryTaxonInfo, to = to, prefix = "ott_id", timeout = timeout
        ))
    }else if( from %in% c("ncbi", "gbif", "worms", "if", "irmng") ){
        x <- paste0(from, ":", x)
        y <- unlist(bplapply(
            x, .queryTaxonInfo, to = to, prefix = "source_id", timeout = timeout
        ))
    }else{
        stop("'from' is not valid.", call. = FALSE)
    }
    # Strip source prefix
    y <- gsub("^.+:", "", y)
    return(y)
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

#' @importFrom jsonlite toJSON
#' @importFrom httr2 req_body_raw req_headers
.queryTNRS <- function(
    x, to, timeout, batch.size = 5000, workers = NULL, factor = 3){
    
    ranges <- .get_batches(x, batch.size, workers, factor)
    
    out.list <- bplapply(
       ranges, function(i) .subqueryTNRS(x[i[1]:i[2]], to, timeout)
    )
    
    out <- do.call(rbind, out.list)
    return(out)
}


.subqueryTNRS <- function(x){
    
    url <- "https://api.opentreeoflife.org/v3/tnrs/match_names"

    body <- list(names = x) |> toJSON(auto_unbox = FALSE)

    resp <- request(url) |>
        req_method("POST") |>
        req_body_raw(charToRaw(body)) |>
        req_headers(`Content-Type` = "application/json") |>
        req_timeout(timeout) |>
        req_perform()
    
    resp <- resp_body_json(resp)

    res <- resp$results
    
    if( to == "ott" ){
        # A
        y <- vapply(
            res,
            function(x){
                if( length(x$matches) == 0L ){
                    NA_integer_
                }else{
                    x$matches[[1]]$taxon$ott_id
                }
            },
            integer(1L)
        )
    }else{
        # B
        y <- vapply(
            res,
            function(x){
                if( length(x$matches) == 0L ){
                    return(NA_character_)
                }
                tax_sources <- unlist(x$matches[[1]]$taxon$tax_sources)
                y <- tax_sources[grepl(to, tax_sources, fixed = TRUE)]
                    
                if (length(y) == 0L) NA_character_ else y
            },
            character(1L)
        )
    }
    return(y)
}
