library(reticulate)

# Use the Python environment
use_python("/Library/Frameworks/Python.framework/Versions/3.9/bin/python3")

# Import python libraries
ssl <- import("ssl")

if( Sys.info()[["sysname"]] == "Darwin" ){
    ssl$`_create_default_https_context` <- ssl$`_create_unverified_context`
}

rdflib <- import("rdflib")

#' @importFrom dplyr
.querySPARQL <- function(module, graph = rdflib$Graph()){
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("\t\t\tuniref:", module, collapse = "\n")
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
    query <- paste0(query_part1, uniref.ids, query_part2)
    # Execute query
    qres <- graph$query(query)
    # Convert name bindings to vector
    tax.vec <- vapply(qres$bindings, function(binding)
        binding[["name"]], character(1))
    return(tax.vec)
}
