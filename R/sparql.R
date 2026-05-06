
#' @importFrom BiocParallel bplapply
#' @importFrom data.table rbindlist
.querySPARQL <- function(from, to, endpoint, init, uniref.identity, timeout,
    batch.size = 25000, workers = NULL, factor = 3){
    
    preface <- paste0("
        PREFIX chebislash: <http://purl.obolibrary.org/obo/chebi/>
        PREFIX enzyme: <http://purl.uniprot.org/enzyme/>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
        PREFIX protein: <http://purl.uniprot.org/uniprot/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        PREFIX up: <http://purl.uniprot.org/core/>
        
        SELECT DISTINCT ?", paste(c(from, to), collapse = " ?"), "
        WHERE {
        "
    )
    # Build triple based on from and to
    triple <- .triple_table(from, to, uniref.identity)
    # Find starts and ends of init ranges
    ranges <- .get_batches(init, batch.size, workers, factor)
    # Query SPARQL for each range
    out.list <- bplapply(ranges, function(i){
        x <- init[i[1]:i[2]]
        query <- .composeSPARQL(from, to, preface, triple, endpoint, x)
        resp <- .sendSPARQL(query, endpoint, timeout)
    })
    # Bind output linkmaps
    out <- rbindlist(out.list)
    return(out)
}


#' @importFrom data.table fread
#' @importFrom httr2 request req_headers req_timeout req_perform resp_body_string req_body_form
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
    linkmap <- fread(text = resp_body_string(resp), header = TRUE)
    return(linkmap)
}


.composeSPARQL <- function(from, to, preface, clause, endpoint, init){
    # If init is defined
    if( !is.null(init) && length(init) != 0L ){
        # Pass values to sciName for taxname ids
        spec.from <- ifelse(from == "taxname", "sciName", from)
        # Limit query with values
        values <- paste0(
            "VALUES ?", spec.from, " {\n",
                paste(.iri_table(init, from, to, endpoint), collapse = " "), "
            }\n")
        # Pre-append values to clause
        clause <- paste0(values, clause)
    }
    # Build query
    query <- paste0(preface, clause, "\n}")
    return(query)
}


.endpoint_table <- function(endpoint){
  
    endpoint <- switch(endpoint,
        Rhea = "https://sparql.rhea-db.org/",
        UniProt = "https://sparql.uniprot.org/"
    )
    
    return(endpoint)
}


.triple_table <- function(from, to, uniref.identity){
    
    spterms <- c(from, to)
    
    internals <- c(
        "uniref", "uniprotkb", "taxname", "taxid", "enzyme", "rhea", "chebi",
        "inchikey", "inchi", "smiles"
    )
    
    ext <- setdiff(spterms, internals)
    
    which_uniref <- ifelse(to == "uniref", paste0("
        ?uniref up:identity ", uniref.identity, ".
    "), "")
    
    uniref2taxid <- paste0(
        which_uniref, "
        # Bind UniRef ids to cluster members
        ?uniref up:member ?member.
        # Bind cluster members to taxa
        ?member up:organism ?taxid.
    ")
    
    taxname2taxid <- "
        # Bind taxa to scientific names and ranks
        ?taxid up:scientificName ?sciName; up:rank ?rank.
    
        # Process name to prefix__taxon format
        BIND(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
        BIND(concat(?prefix, '__', ?sciName) as ?taxname)
    "
    
    uniref2taxname <- paste0(uniref2taxid, taxname2taxid)
    
    uniref2uniprotkb <- paste0(
        which_uniref, "
        ?uniprotkb up:representativeFor ?uniref.
    ")
    
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
        ?compound rh:chebi|rh:reactivePart/rh:chebi|rh:underlyingChebi ?chebi.
    "
    
    inchi2chebi <- "
        ?chebi chebislash:inchi ?inchi.
    "
    
    inchikey2chebi <- "
        ?chebi chebislash:inchikey ?inchikey.  
    "
    
    smiles2chebi <- "
        ?chebi chebislash:smiles ?smiles.
    "
    
    rhea2inchi <- paste0(rhea2chebi, inchi2chebi)
    rhea2inchikey <- paste0(rhea2chebi, inchikey2chebi)
    smiles2rhea <- paste0(rhea2chebi, smiles2chebi)
    
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
        FILTER(contains(str(?", ext, "), '", ext, "'))
    ")
    
    chebi2external <- paste0("
        ?chebi oboInOwl:hasDbXref ?", ext, ".
        FILTER(contains(?", ext, ", '", ext, ":'))
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
        inchi2chebi = inchi2chebi,
        inchikey2chebi = inchikey2chebi,
        smiles2chebi = smiles2chebi,
        rhea2inchi = rhea2inchi,
        rhea2inchikey = rhea2inchikey,
        smiles2rhea = smiles2rhea,
        "external"
    )
    
    if( triple == "external" ){
      
        internal <- setdiff(spterms, ext)
        
        triple <- switch(
            internal,
            uniprotkb = uniprotkb2external,
            uniref = uniref2external,
            rhea = rhea2external,
            chebi = chebi2external
        )
    }
    
    return(triple)
}


.iri_table <- function(x, from, to, endpoint) {
    
    iri <- switch(
        from,
        chebi = "obo",
        uniprotkb = "protein",
        uniref = "uniref",
        taxid = "taxon",
        smiles = "",
        inchi = "",
        inchikey = "",
        taxname = "",
        enzyme = "enzyme",
        rhea = "rh",
        "external"
    )
    
    is_external <- iri == "external"
    
    if( is_external && to == "chebi" ){
        x <- paste0("'", from, ":", x, "'")
    }else if( is_external ){
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
            FILTER(contains(str(?entry), '", ext, "'))
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
