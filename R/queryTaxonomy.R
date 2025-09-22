#' Query taxonomy from UniProt using SPARQL
#' 
#' @name queryTaxonomy
NULL

#' @export
#' @rdname queryTaxonomy
#' @importFrom BiocParallel bplapply
setMethod("queryTaxonomy", signature = c(map = "list"),
    function(map, remove.empty = TRUE, verbose = TRUE){
        if( !is.logical(remove.empty) ){
            stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
        }
        if( !is.logical(verbose) ){
            stop("'verbose' should be TRUE or FALSE.", call. = FALSE)
        }
        # Query taxonomy from UniProt
        sig.list <- bplapply(map, .querySPARQL)
        if( verbose ){
            message(length(unlist(sig.list)), " taxa queried from UniProt.")
        }
        # Remove empty modules
        if( remove.empty ){
            sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
        }
        return(sig.list)
    }
)


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