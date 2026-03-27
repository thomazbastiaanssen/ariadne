
#' @importFrom BiocParallel bplapply
.querySPARQL <- function(from, to, endpoint, init, timeout,
    batch.size = 25000, workers = NULL, factor = 3){
    
    ranges <- .get_batches(init, batch.size, workers, factor)
    
    out.list <- bplapply(ranges, function(i){
        x <- init[i[1]:i[2]]
        query <- .composeSPARQL(from, to, endpoint, x)
        resp <- .sendSPARQL(query, endpoint, timeout)
    })
    
    out <- do.call(rbind, out.list)
    return(out)
}


#' @importFrom utils read.csv
#' @importFrom httr2 request req_headers req_timeout req_perform
#'   resp_body_string req_body_form
.sendSPARQL <- function(query, endpoint, timeout) {
    # Get base url from endpoint table
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


.composeSPARQL <- function(from, to, endpoint, init) {
    
    preface <- "
        PREFIX enzyme: <http://purl.uniprot.org/enzyme/>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX protein: <http://purl.uniprot.org/uniprot/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        PREFIX up: <http://purl.uniprot.org/core/>
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
                paste0(.iri_table(init, from, endpoint), collapse = " "), "\n",
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
        Rhea = "https://sparql.rhea-db.org/",
        UniProt = "https://sparql.uniprot.org/"
    )
    
    return(endpoint)
}


.triple_table <- function(from, to){
    
    spterms <- c(from, to)
    
    internals <- c(
        "uniref", "uniprotkb", "taxname", "taxid", "enzyme", "rhea", "chebi"
    )
    
    ext <- setdiff(spterms, internals)
    
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
        ?uniprotkb (up:enzyme|up:domain/up:enzyme|up:component/up:enzyme) ?enzyme.
    "
    
    uniref2enzyme <- paste0(uniref2uniprotkb, uniprotkb2enzyme)
    
    uniprotkb2external <- paste0("
        ?uniprotkb a up:Protein.
        ?uniprotkb rdfs:seeAlso ?", ext, ".
        ?", ext, " up:database <http://purl.uniprot.org/database/", ext, ">.
    ")
    
    uniref2external <- paste0(uniref2uniprotkb, uniprotkb2external)
    
    uniprotkb2rhea <- "
        ?uniprotkb up:annotation/up:catalyticActivity/up:catalyzedReaction ?rhea.
    "
    
    uniref2rhea <- paste0(uniref2uniprotkb, uniprotkb2rhea)
    
    rhea2chebi <- "
        ?rhea rdfs:subClassOf rh:Reaction.
        ?rhea rh:side/rh:contains/rh:compound ?compound.
        # chebi is small molecule, reactive part of macromolecule or polymer
        ?compound (rh:chebi|(rh:reactivePart/rh:chebi)|rh:underlyingChebi) ?chebi.
    "
    
    rhea2enzyme <- "
        ?rhea rdfs:subClassOf rh:Reaction.
        ?rhea rh:ec ?enzyme.
    "
    
    rhea2external <- paste0("
        ?master rdfs:subClassOf rh:Reaction.
        {
            ?master rdfs:seeAlso ?", ext, ".
            BIND(?master as ?rhea)
        }
        UNION
        {
            ?master rh:directionalReaction|rh:bidirectionalReaction ?alter.
            ?alter rdfs:seeAlso ?", ext, ".
            BIND(?alter as ?rhea)
        }
        FILTER(CONTAINS(str(?", ext, "), '", ext, "'))
    ")
    
    key <- paste(rev(sort(spterms)), collapse = "2")
    
    triple <- switch(
        key,
        uniref2taxid = uniref2taxid,
        taxname2taxid = taxname2taxid,
        uniref2taxname = uniref2taxname,
        uniref2uniprotkb = uniref2uniprotkb,
        uniref2enzyme = uniref2enzyme,
        uniprotkb2enzyme = uniprotkb2enzyme,
        uniprotkb2rhea = uniprotkb2rhea,
        uniref2rhea = uniref2rhea,
        rhea2chebi = rhea2chebi,
        rhea2enzyme = rhea2enzyme,
        "external"
    )
    
    if( triple == "external" ){
      
        internal <- setdiff(spterms, ext)
        
        triple <- switch(
            internal,
            uniprotkb = uniprotkb2external,
            uniref = uniref2external,
            rhea = rhea2external
        )
    }
    
    return(triple)
}


.iri_table <- function(x, from, endpoint) {
    
    iri <- switch(
        from,
        chebi = "obo",
        uniprotkb = "protein",
        uniref = "uniref",
        taxid = "taxon",
        taxname = "",
        enzyme = "enzyme",
        rhea = "rh",
        "external"
    )
    
    if( iri == "external" ){
        iri <- .get_external_iri(from, endpoint)
        x <- paste0("<", iri, x, ">")
    }else if( iri == "" ){
        x <- paste0("'", x, "'")
    }else{
        x <- paste0(iri, ":", x)
    }
    
    return(x)
}


.get_external_iri <- function(ext, endpoint, limit = 10) {
    
    uniprot_query <- paste0("
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?entry
        WHERE {
            ?protein a up:Protein.
            ?protein rdfs:seeAlso ?entry.
            ?entry up:database <http://purl.uniprot.org/database/", ext, ">.
    ")
    
    rhea_query <- paste0("
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?entry
        WHERE {
            ?rhea rdfs:subClassOf rh:Reaction.
            ?rhea rdfs:seeAlso ?entry.
            FILTER CONTAINS(str(?entry), '", ext, "')
    ")
    
    query <- switch(endpoint, Rhea = rhea_query, UniProt = uniprot_query)
    query <- paste0(query, "} LIMIT ", limit)
    
    iri_list <- .sendSPARQL(query, endpoint, 1e6)
    iri_list <- unlist(iri_list, use.names = FALSE)
    iri_list <- gsub("([^/]+)$", "", iri_list)
    iri <- unique(iri_list)
    
    if( length(iri) > 1L ){
        stop("Multiple prefixes were found in the selected database",
            call. = FALSE)
    }
    
    return(iri)
}


#' @importFrom BiocParallel bpworkers
.get_batches <- function(x, batch.size, workers, factor){
    
    if( is.null(x) ){
        return(list(c(1L, 1L)))
    }
    
    xlen <- length(x)
    
    if( is.null(workers) ){
        workers <- bpworkers()
    }
    
    batch.num <- min(ceiling(xlen / batch.size), factor * workers)
    adapted.size <- ceiling(xlen / batch.num)
    
    if( adapted.size > batch.size ){
        stop("Query limit was reached (", adapted.size, " > ", batch.size, ").",
            " Increase 'factor', 'batch.size' or 'workers' and try again.",
            call. = FALSE)
    }
    
    ranges <- lapply(seq_len(batch.num), function(i) {
        start <- (i - 1) * adapted.size + 1
        end <- min(i * adapted.size, xlen)
        c(start, end)
    })
    
    return(ranges)
}
