library(reticulate)

# Use the Python environment
use_python("/Library/Frameworks/Python.framework/Versions/3.9/bin/python3")

# Import python libraries
ssl <- import("ssl")

if( Sys.info()[["sysname"]] == "Darwin" ){
    ssl$`_create_default_https_context` <- ssl$`_create_unverified_context`
}

rdflib <- import("rdflib")

.querySPARQL <- function(uniref.ids, graph = rdflib$Graph()){
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("uniref:", uniref.ids, collapse = " ")
    # Define first part of query
    query_part1 <- "
        PREFIX uniprot: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        # Outer query to get final names
        SELECT ?unirefOut ?name
        WHERE {
            SERVICE <https://sparql.uniprot.org/> {
            {   # Inner query to get distinct taxa
                SELECT DISTINCT ?unirefIn ?taxId
                WHERE {
                    # List UniRef90 ids
                    VALUES ?unirefIn {
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
            # Strip prefix from uniref ids
            BIND(strafter(str(?unirefIn), 'uniref/') as ?unirefOut)
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
        binding[["name"]], character(1))
    # Convert uniref bindings to vector
    uniref.vec <- vapply(qres$bindings, function(binding)
        binding[["unirefOut"]], character(1))
    # Combine into linkMap
    uniref.tax <- data.frame(uniref = uniref.vec, tax = tax.vec)
    return(uniref.tax)
}
