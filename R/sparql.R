library(reticulate)

# Use the Python environment
use_python("/Library/Frameworks/Python.framework/Versions/3.9/bin/python3")

# Import rdflib
ssl <- import("ssl")
rdflib <- import("rdflib")

ssl$`_create_default_https_context` <- ssl$`_create_unverified_context`

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
            ?unirefIn uniprot:member ?member .
            ?member uniprot:organism ?taxIn .
            ?taxIn uniprot:scientificName ?sciName .
            ?taxIn uniprot:rank ?rank .
  
            BIND(STRAFTER(STR(?unirefIn), '_') AS ?unirefId)
            BIND(STRAFTER(STR(?taxIn), 'taxonomy/') AS ?taxId)

            BIND(LCASE(SUBSTR(STRAFTER(STR(?rank), 'Rank_'), 1, 1)) AS ?prefix)
            BIND(CONCAT(?prefix, '__', ?sciName) AS ?name)
            }
        }
        "
    
    query <- paste0(query_part1, uniref.ids, query_part3)
    
    # Execute the query
    qres <- graph$query(query)
    
    # Convert the bindings to a data frame
    res.df <- bind_rows(lapply(qres$bindings, function(binding) {
        data.frame(
            unirefId = binding[["unirefId"]],
            taxId = binding[["taxId"]],
            name = binding[["name"]]
        )
    }))
    
    return(res.df)
}