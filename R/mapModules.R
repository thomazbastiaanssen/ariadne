#' Map modules to taxa
#' 
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
#' @param type \code{Character scalar}. Specifies the type of modules to expect
#'   as input. It can be one of \code{c("oto", "andor")}.
#'   (Default: \code{"oto"}).
#' 
#' @param mode \code{Character scalar}. Specifies the type of modules to return
#'   as output. It can be one of \code{c("uniref", "taxonomy")}.
#'   (Default: \code{"uniref"}).
#'
#' @param remove.empty \code{Logical scalar}. Should modules with no matching
#'   taxa be removed. (Default: \code{TRUE}).
#' 
#' @param verbose \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#' 
#' @details
#' The input modules of \code{mapModules} can be either one-to-one mapped or
#' and/or relationships, depending on whether \code{type} is \code{"oto"} or
#' \code{"andor"}. The output modules of \code{mapModules} can be either uniref
#' or taxonomies, depending on whether \code{mode} is \code{"uniref"} or
#' \code{"taxonomy"}.
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
#' map <- importMapping("ChocoPhlAn", from = "ko", to = "uniref90")
#' 
#' # Map modules to UniRef90
#' uniref.modules <- mapModules(gbm[seq(3)], map)
#' 
#' # Map modules to taxa
#' tax.modules <- mapModules(gbm[seq(3)], map, mode = "taxonomy")
#' 
#' @name mapModules
NULL

#' @rdname mapModules
#' @export
#' @importFrom BiocParallel bplapply
setMethod("mapModules", signature = c(modules = "list"),
    function(
        modules, map, type = "oto", mode = "uniref", remove.empty = TRUE,
        verbose = TRUE
    ){
        # Check arguments
        if( !is.vector(modules) ){
            stop("'modules' must be a character vector or list of character ",
                "vectors, where each vector corresponds to a module.",
                call. = FALSE)
        }
        if( !is.vector(map) ){
            stop("'map' must be a character vector or list of character ",
                "vectors, where each vector corresponds to a mapping.",
                call. = FALSE)
        }
        if( !mode %in% c("uniref", "taxonomy") ){
            stop("'mode' must be either uniref or taxonomy.", call. = FALSE)
        }
        if( !is.logical(remove.empty) ){
            stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
        }
        # Keep only relevant bindings
        keep <- names(map) %in% unique(unlist(modules, use.names = FALSE))
        map <- map[keep]
        # Check matched uniref ids
        map.size <- length(unlist(map))
        if( map.size == 0){
            stop("'map' did not match any element in 'modules'.", call. = FALSE)
        }
        if( verbose ){
            message("Mapping ", length(modules), " modules to ",
                map.size, " values.")
        }
        if( mode == "taxonomy" ){
            # Query taxonomy from uniref
            map <- bplapply(map, .querySPARQL)
        }
        # Select mapper based on module type
        mapper <- switch(type,
            oto = .oto_mapping,
            andor = .andor_mapping
        )
        # Map modules and store in modules list
        sig.list <- mapper(modules, map)
        # Remove empty modules
        if( remove.empty ){
            sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
        }
        return(sig.list)
    }
)

# Query taxonomies from uniref ids from UniProt using SPARQL
.querySPARQL <- function(module, graph = rdflib$Graph()){
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("uniref:", module, collapse = " ")
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

# Perform one-to-one mapping
.oto_mapping <- function(modules, values){
    # Store taxa in modules list
    sig.list <- bplapply(modules, function(module){
        keep <- names(values) %in% module
        module <- unname(unlist(values[keep]))
        return(module)
    })
    return(sig.list)
}

# Perform and/or mapping (reaction pathway modules)
.andor_mapping <- function(modules, values){
    # Find functions for each taxon
    linkmap <- as.linkmap(values)
    tax <- split(linkmap$x, linkmap$y)
    # Store taxa in modules list
    sig.list <- bplapply(modules, function(module) {
        members <- vapply(tax, function(tax.item) {
            all(vapply(module, function(comp) any(comp %in% tax.item), logical(1L)))
        }, logical(1L))
        names(tax)[members]
    })
    return(sig.list)
}
