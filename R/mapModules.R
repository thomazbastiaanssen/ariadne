#' Map modules to taxa
#' @name mapModules
#' @rdname mapModules
#'
#' @description
#' \code{mapModules} returns a list of modules containing the taxa that are
#' members of each module. Taxa are derived from uniref ids by querying UniProt
#' SPARQL. Membership is based on whether a taxon meets the criteria specified
#' by a module.
#'
#' @param modules \code{Named character list}. Named list of vectors, where each
#'   vector is a module and its elements are the module members.
#'
#' @param map \code{Named character list}. Named list of vectors, where each
#'   vector is a mapping key and its elements are the mapped values.
#'
#' @param mode \code{Character scalar}. Specifies the type of modules to expect
#'   as input. It can be one of \code{c("single", "andor")}.
#'   (Default: \code{"single"}).
#'
#' @param uniprot \code{Logical scalar}. Should a SPARQL query to UniProt be
#'   made to map UniRef ids to taxa (Default: \code{TRUE}).
#'
#' @param remove.empty \code{Logical scalar}. Should modules with no matching
#'   taxa be removed. (Default: \code{TRUE}).
#'
#' @param verbose \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#'
#' @details
#' The input modules of \code{mapModules} can be either single elements or
#' and/or relationships, depending on whether \code{mode} is \code{"single"} or
#' \code{"andor"}.
#'
#' @return
#' \code{mapModules} returns a named list of vectors, where each vector is a
#' module and its elements are the members of that module.
#'
#' @examples
#' # Import GBM
#' gbm <- importModules("GBM")
#'
#' # Import ko-to-uniref90 mapping
#' map <- importMapping(ChocoPhlAn, ko ~ uniref90, dry_run = FALSE)
#'
#'
#' x <- c(gbm, map)
#' # Map modules to UniRef90
#' modules <- mapModules(
#'     x,
#'     mode = "andor",
#'     uniprot = TRUE
#' )
NULL


#' @importFrom BiocParallel bplapply
S7::method(mapModules, MultiFactor) <- function(
        x, mapping = module ~ uniref90, verbose = TRUE, dry_run = TRUE
        ) {

    # Check arguments
    if( !is.logical(uniprot) ){
        stop("'uniprot' should be TRUE or FALSE.", call. = FALSE)
    }
    if( !is.logical(verbose) ){
        stop("'verbose' should be TRUE or FALSE.", call. = FALSE)
    }
    linkmap <- weave(x, mapping)
    query   <- as.character(unique(linkmap[[2L]]))

    if( uniprot ){
        # Query taxonomy from UniProt
        query <- bplapply(query, .querySPARQL)
    }
    if( uniprot && verbose ){
        message(length(unlist(query)), " taxa were queried from UniProt.")
    }
    # Select mapping method based on module type
    map.method <- switch(mode,
        single = .single_mapping,
        andor = .andor_mapping
    )
    # Map modules and store in modules list
    sig.list <- map.method(modules, map)
    # Remove empty modules
    if( remove.empty ){
        sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
    }
    if( verbose ){
        message(length(unlist(sig.list)), " items were mapped to ",
            length(sig.list), " modules.")
    }
    return(sig.list)
}

# Perform one-to-one mapping
.single_mapping <- function(x, y){
    # Filter y keys that match x values
    y <- y[names(y) %in% unique(unlist(x, use.names = FALSE))]
    if( length(y) == 0 ){
        stop("Items in 'x' did not match any item in 'y'.", call. = FALSE)
    }
    # Store taxa in modules list
    z <- bplapply(x, function(values)
        unlist(y[names(y) %in% values], use.names = FALSE))
    return(z)
}

#' @importFrom MultiFactor as.LinkMap
# Perform and/or mapping (reaction pathway modules)
.andor_mapping <- function(x, y){
    # Filter y keys that match x values
    y <- y[names(y) %in% unique(unlist(x, use.names = FALSE))]
    # Find functions for each taxon
    linkmap <- as.LinkMap(y)
    tax <- split(linkmap$x, linkmap$y)
    # Store taxa in modules list
    z <- lapply(x, function(module) {
        members <- vapply(tax, function(tax.item) {
            all(vapply(module, function(comp) any(comp %in% tax.item), logical(1L)))
        }, logical(1L))
        names(tax)[members]
    })
    return(z)
}

# Query taxonomies based on uniref ids from UniProt using SPARQL
.querySPARQL <- function(x, graph = rdflib$Graph()){
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("uniref:", x, collapse = " ")
    # Define first part of query
    query_part1 <- "
        PREFIX uniprot: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        # Outer query to get final names
        SELECT ?name
        WHERE {
            SERVICE <https://sparql.uniprot.org/> {
            {   # Inner query to get distinct taxa
                SELECT DISTINCT ?taxId
                WHERE {
                    # List UniRef90 ids
                    VALUES ?unirefId {
        "
    # Define second part of query
    query_part2 <- "
                    }
                    # Bind UniRef90 ids to cluster members
                    ?unirefId uniprot:member ?member.
                    # Bind cluster members to taxa
                    ?member uniprot:organism ?taxId.
                }
            }
            # Bind taxa to scientific names and ranks
            ?taxId uniprot:scientificName ?sciName; uniprot:rank ?rank.
            # Process name to prefix__taxon format
            BIND(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
            BIND(concat(?prefix, '__', ?sciName) as ?name)
            }
        }
        "
    # Build query
    query <- paste0(query_part1, "\t\t\t", uniref.ids, query_part2)
    # Execute query
    qres <- graph$query(query)
    # Convert name bindings to vector
    tax.vec <- vapply(qres$bindings, function(binding)
        binding[["name"]], character(1L))
    return(tax.vec)
}
