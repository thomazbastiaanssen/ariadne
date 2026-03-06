
#' @importFrom BiocParallel bplapply
.querySPARQL <- function(from, to, endpoint, init, timeout,
    batch.size = 25000, workers = NULL, factor = 3){
    
    ranges <- .get_batches(init, batch.size, workers, factor)
    
    out.list <- bplapply(ranges, function(i){
        x <- init[i[1]:i[2]]
        query <- .composeSPARQL(from, to, x)
        resp <- .sendSPARQL(query, endpoint, timeout)
    })
    
    out <- do.call(rbind, out.list)
    colnames(out) <- c(from, to)
    return(out)
}


#' @importFrom utils URLencode read.csv
#' @importFrom httr2 request req_headers req_timeout req_perform
#'    resp_check_status resp_body_string req_body_form
.sendSPARQL <- function(query, endpoint, timeout) {
    # Get base URL from endpoint table
    endpoint <- .endpoint_table(endpoint)
    # Build request with Accept header for CSV format
    req <- request(endpoint) |>
        req_method("POST") |>
        req_body_form(query = query) |>
        req_headers(Accept = "text/csv") |>
        req_timeout(timeout)
    # Get response
    resp <- req_perform(req)
    # Parse CSV content into data frame
    linkmap <- read.csv(
        text = resp_body_string(resp),
        stringsAsFactors = FALSE
    )
    return(linkmap)
}


.composeSPARQL <- function(from, to, init) {
    
    preface <- "
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    "

    clause <- paste0("
        SELECT DISTINCT ?", paste0(c(from, to), collapse = " ?"), "
        WHERE {
        "
    )
    
    if( !is.null(init) && length(init) != 0L ){
        
        spec.from <- ifelse(from == "taxname", "sciName", from)
        
        clause <- paste0(
            clause,
            "VALUES ?", spec.from, " {\n",
                paste0(.iri_table(init, from), collapse = " "), "\n",
            "}\n"
        )
    }
    
    triple <- .triple_table(from, to)
    
    clause <- paste0(clause, triple, "\n}")
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
    
    uniprotkb2enzyme <- "
        ?uniprotkb ( up:enzyme | up:domain/up:enzyme | up:component/up:enzyme ) ?enzyme.
    "
    
    uniref2enzyme <- paste0(uniref2uniprotkb, uniprotkb2enzyme)
    
    ext <- spterms[!spterms %in% c("uniref", "uniprotkb", "taxname", "taxid", "enzyme")]

    uniprotkb2external <- paste0("
        ?uniprotkb a up:Protein.
        ?uniprotkb rdfs:seeAlso ?", ext, ".
        ?", ext, " up:database <http://purl.uniprot.org/database/", ext, ">.
    ")
    
    uniref2external <- paste0(uniref2uniprotkb, uniprotkb2external)
    
    key <- paste(rev(sort(spterms)), collapse = "2")
    
    triple <- switch(
        key,
        uniref2taxid = uniref2taxid,
        taxname2taxid = taxname2taxid,
        uniref2taxname = uniref2taxname,
        uniref2uniprotkb = uniref2uniprotkb,
        uniref2enzyme = uniref2enzyme,
        uniprotkb2enzyme = uniprotkb2enzyme,
        "external"
    )
    
    if( triple == "external" && "uniprotkb" %in% spterms ){
        triple <- uniprotkb2external
    }else if( triple == "external" && "uniref" %in% spterms ){
        triple <- uniref2external
    }
    
    return(triple)
}


.iri_table <- function(x, from) {
    
    iri <- switch(
        from,
        uniprotkb = "up",
        uniref = "uniref",
        taxid = "taxon",
        taxname = "",
        enzyme = "enzyme",
        "external"
    )

    if( iri == "external" ){
        iri <- .get_external_iri(from)
        x <- paste0("<", iri, x, ">")
    }else if( iri == "" ){
        x <- paste0("'", x, "'")
    }else{
        x <- paste0(iri, ":", x)
    }
    
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
    
    iri_list <- .sendSPARQL(query, "UniProt", 1e6)
    iri_list <- unlist(iri_list, use.names = FALSE)
    iri_list <- gsub("([^/]+)$", "", iri_list)
    iri <- unique(iri_list)
    
    if( length(iri) > 1L ){
        stop("Multiple prefixes were found in the selected database",
            call. = FALSE)
    }
    
    return(iri)
}

