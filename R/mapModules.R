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
#' @param remove.empty \code{Logical scalar}. Should modules with no matching
#'   taxa be removed. (Default: \code{TRUE}).
#' 
#' @param message \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#' 
#' @details
#' Additional details
#' 
#' @return
#' \code{mapModules} returns a named list of vectors, where each vector is a
#' module and its elements are the taxa members of that module.
#'
#' @examples
#' # Import GBM
#' gbm <- importModules("GBM")
#' 
#' # Import ko-to-uniref90 mapping
#' map <- importMapping("ChocoPhlAn", from = "ko", to = "uniref90")
#' 
#' # Map modules to taxa
#' sigs <- mapModules(gbm[seq(3)], map)
#' 
#' @name mapModules
NULL

#' @rdname mapModules
#' @export
#' @importFrom BiocParallel bplapply
setMethod("mapModules", signature = c(modules = "list"),
    function(modules, map, remove.empty = TRUE, message = TRUE){
        # Check arguments
        if( !is.vector(modules) ){
          stop("'modules' must be a character vector or list of character ",
              "vectors, where each vector corresponds to a module.", call. = FALSE)
        }
        if( !is.vector(map) ){
          stop("'map' must be a character vector or list of character ",
              "vectors, where each vector corresponds to a mapping.", call. = FALSE)
        }
        if( !is.logical(remove.empty) ){
            stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
        }
        # Keep only relevant bindings
        keep <- names(map) %in% unique(unlist(modules, use.names = FALSE))
        map <- map[keep]
        # Check matched uniref ids
        uniref.length <- length(unlist(map))
        if( uniref.length == 0){
            stop("'map' did not match any element in 'modules'.", call. = FALSE)
        }
        if( message ){
            message("Mapping ", length(modules), " modules to ",
                uniref.length, " uniref ids.")
        }
        # Query taxonomy from UniRef90
        tax.list <- bplapply(map, .querySPARQL)
        # Find functions for each taxon
        tax.linkmap <- as.linkmap(tax.list)
        tax <- split(tax.linkmap$x, tax.linkmap$y)
        if( message ){
            message(length(tax), " tax ids were retrieved from UniProt.")
        }
        # Store taxa in modules list
        sig.list <- bplapply(modules, function(module) {
            members <- vapply(tax, function(tax.item) {
                all(vapply(module, function(mod) any(mod %in% tax.item), logical(1)))
            }, logical(1))
            names(tax)[members]
        })
        if( remove.empty ){
            # Remove empty modules 
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
            bind(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
            bind(concat(?prefix, '__', ?sciName) as ?name)
            }
        }
        "
    # Build query
    query <- paste0(query_part1, "\t\t\t", uniref.ids, query_part2)
    # Execute query
    qres <- graph$query(query)
    # Convert name bindings to vector
    tax.vec <- vapply(qres$bindings, function(binding)
        binding[["name"]], character(1))
    return(tax.vec)
}
