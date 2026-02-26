#' @title SPARQLmap
#' @name SPARQLmap
#' @param x `Character vector`. Optional. Feature IDs on which the query should
#'     be constrained.
#' @param x.type `Character scalar`. Feature type of `x`. (Default: `"uniref"`)
#' @param .by either a `formula`, or a `Character vector` of length 2,  with the
#'      names of the desired combination of feature types.
#' @param endpoint The SPARQL endpoint (a URL)
#' @returns a data.frame of desired query.
#' @examples
#'
#'
#'x <- paste0(
#'    "UniRef90_",
#'    c("A0A010PZR5", "A0A010PZT4", "A0A010PZU0", "A0A010PZV8",
#'    "A0A010PZW7", "A0A010Q006", "A0A010Q047", "A0A010Q0B1")
#'    )
#'
#'
#' SPARQLmap(x, .by = uniref ~ ec,  endpoint = "uniprot")
#'
#' SPARQLmap(x, .by = uniref ~ species,  endpoint = "uniprot")
#'
.querySPARQL <- function(x, from, to, endpoint, timeout = 1e6){
    
    query <- .composeSPARQL(x, from, to)
    
    out <- .sendSPARQL(query, endpoint, timeout)
    
    colnames(out) <- c(from, to)
    
    return(out)
}


#' @importFrom utils URLencode read.csv
#' @importFrom httr2 request req_headers req_timeout req_perform
#'    resp_check_status resp_body_string
#' @noRd
.sendSPARQL <- function(query, endpoint, timeout = 1e6) {
    # Get base URL from endpoint table
    endpoint <- .endpoint_table(endpoint)
    # URL-encode only the query string (the SPARQL query)
    query <- URLencode(query, reserved = TRUE)
    # Construct full URL by appending query parameter
    query <- paste0(endpoint, "?query=", query)
    # Build request with Accept header for CSV format
    req <- request(query) |>
        req_headers(Accept = "text/csv") |>
        req_timeout(timeout)
    # Get response
    resp <- req_perform(req)
    # Check for HTTP errors
    resp_check_status(resp)
    # Parse CSV content into data frame
    linkmap <- read.csv(
        text = resp_body_string(resp),
        stringsAsFactors = FALSE
    )
    return(linkmap)
}


#' @importFrom utils URLencode
.composeSPARQL <- function(x, from, to) {
    
    if( length(x) == 0L ){
        stop("No binding found", call. = FALSE)
    }
    
    preface <- "
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    "
    
    triple <- .triple_table(from, to)

    clause <- paste0("
        SELECT DISTINCT ?", paste0(c(from, to), collapse = " ?"), "
        WHERE {
            VALUES ?", from, " {\n",
                paste0(.iri_table(x, from), collapse = " "), "\n",
            "}",
            triple,
        "}"
    )
    # Build query
    query <- paste0(preface, clause)
    return(query)
}


.endpoint_table <- function(endpoint){
  
    endpoint <- switch(endpoint,
        UniProt = "https://sparql.uniprot.org/"
    )
    
    return(endpoint)
}


.triple_table <- function(from, to){
    
    spterms <- c(from, to)
    
    uniref2taxid <- "
        # Bind UniRef ids to cluster members
        ?uniref up:member ?member.
        # Bind cluster members to taxa
        ?member up:organism ?taxid.
    "
    
    taxname2taxid <- "
        # Bind taxa to scientific names and ranks
        ?taxid up:scientificName ?sciName; up:rank ?rank.
    
        # Process name to prefix__taxon format
        BIND(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
        BIND(concat(?prefix, '__', ?sciName) as ?taxname)
    "
    
    uniref2taxname <- paste0(uniref2taxid, taxname2taxid)
    
    uniref2uniprotkb <- "
        ?uniprotkb up:representativeFor ?uniref.
    "
    
    uniprotkb2ec <- "
        ?uniprotkb ( up:enzyme | up:domain/up:enzyme | up:component/up:enzyme ) ?ec.
    "
    
    uniref2ec <- paste0(uniref2uniprotkb, uniprotkb2ec)
    
    external <- spterms[!spterms %in% c("uniref", "uniprotkb", "taxid", "taxname", "ec")]
    
    uniprotkb2external <- paste0("
        ?uniprotkb a up:Protein.
        ?uniprotkb rdfs:seeAlso ?", external, ".
        ?", external, " up:database <http://purl.uniprot.org/database/", external, ">.
    ")
    
    uniref2external <- paste0(uniref2uniprotkb, uniprotkb2external)
    
    key <- paste(rev(sort(spterms)), collapse = "2")
    
    triple <- switch(
        key,
        uniref2taxid = uniref2taxid,
        taxname2taxid = taxname2taxid,
        uniref2taxname = uniref2taxname,
        uniref2uniprotkb = uniref2uniprotkb,
        uniref2ec = uniref2ec,
        uniprotkb2ec = uniprotkb2ec,
        "external"
    )
    
    if( triple == "external" && "uniprotkb" %in% spterms){
        triple <- uniprotkb2external
    }else if( triple == "external" && "uniref" %in% spterms ){
        triple <- uniref2external
    }
    
    return(triple)
}


.iri_table <- function(x, from) {
    
    iri <- switch(
        from,
        "uniref" = "uniref",
        "taxid" = ,
        "taxname" = ,
        "uniprotkb" = "up",
        "ec" = "enzyme",
        "external"
    )

    if( iri == "external" ){
        iri <- .get_external_iri(from)
        x <- paste0("<", iri, x, ">")
    }else(
        x <- paste0(iri, ":", x)
    )
    
    return(x)
}



.get_external_iri <- function(ext, limit = 10) {
    
    query <- paste0("
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?entry
        WHERE {
            ?protein a up:Protein.
            ?protein rdfs:seeAlso ?entry.
            ?entry up:database <http://purl.uniprot.org/database/", ext, ">.\n",
        "}
        LIMIT ", limit
    )
    
    iri_list <- .sendSPARQL(query, "UniProt")
    iri_list <- unlist(iri_list, use.names = FALSE)
    iri_list <- gsub("([^/]+)$", "", iri_list)
    iri <- unique(iri_list)
    
    if( length(iri) != 1L ){
        stop("Multiple prefixes were found in the selected database",
            call. = FALSE)
    }
    
    return(iri)
}

