library(reticulate)

# Use the Python environment
use_python("/Library/Frameworks/Python.framework/Versions/3.9/bin/python3")

# Import python libraries
ssl <- import("ssl")

if( Sys.info()[["sysname"]] == "Darwin" ){
    ssl$`_create_default_https_context` <- ssl$`_create_unverified_context`
}

rdflib <- import("rdflib")

#' @importFrom dplyr bind_rows
.querySPARQL <- function(module, graph = rdflib$Graph()){
    
    uniref.ids <- paste0("\t\t\tuniref:", module, collapse = "\n")
    
    query_part1 <- "
        PREFIX uniprot: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT ?unirefId ?taxId ?name
        WHERE {
            SERVICE <https://sparql.uniprot.org/> {
                VALUES ?unirefIn {
        "
    
    query_part3 <- "
                }
            ?unirefIn uniprot:member ?member.
            ?member uniprot:organism ?taxIn.
            ?taxIn uniprot:scientificName ?sciName; uniprot:rank ?rank.
  
            bind(strafter(str(?unirefIn), '_') as ?unirefId)
            bind(strafter(str(?taxIn), 'taxonomy/') as ?taxId)

            bind(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
            bind(concat(?prefix, '__', ?sciName) as ?name)
            }
        }
        "
    
    # Build query
    query <- paste0(query_part1, uniref.ids, query_part3)
    
    # Execute query
    qres <- graph$query(query)
    
    # Convert name bindings to vector
    tax.vec <- vapply(qres$bindings, function(binding)
        binding[["name"]], character(1))
    
    return(tax.vec)
}
