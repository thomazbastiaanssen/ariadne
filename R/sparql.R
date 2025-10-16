
sparkmap <- function(x, endpoint = "uniprot", prefixes) {
    known_endpoints <- c("uniprot")
    endpoint <- known_endpoints[charmatch(endpoint, table = known_endpoints)]
    endpoint_url <- switch(endpoint,  "uniprot" = "https://sparql.uniprot.org/")

    query <- .composeSPARQL(x, endpoint_url, prefixes)
    .querySPARQL(query)
}

#' @param endpoint The SPARQL endpoint (a URL)
#' @importFrom utils URLencode
.composeSPARQL <- function(x, endpoint_url, prefixes) {
    #uniref.ids <- paste0("uniref:", x, collapse = " ")


    #TODO translate ariadne IDs to uniref code.

    # Translate ariadne path to SPARQL
    # (paste rows sequentially)


    q <- paste0("
SELECT ?unirefId ?sciName ?member ?protein ?enzyme
WHERE {
    VALUES ?unirefId {
    ", paste0("uniref:", x, collapse = " "),
                "
    } .
?unirefId up_core:member ?member .
?member up_core:organism ?taxId .
?taxId up_core:scientificName ?sciName .
?protein up_core:representativeFor ?unirefId .
?protein ( up_core:enzyme | up_core:domain/up_core:enzyme | up_core:component/up_core:enzyme ) ?enzyme .
}")
    # Build query
    query <- paste0(prefixes, q)
    query <- gsub("\\+", "%2B", utils::URLencode(query, reserved = TRUE))
    query <- paste0(endpoint_url, "?query=", query)
    return(query)
}


#' @param query The SPARQL query (character)
#' @returns a data.frame of desired query.
#' @examples
#' endpoint <- "https://lindas.admin.ch/query"
#' query <- "PREFIX schema: <http://schema.org/>
#'   SELECT * WHERE {
#'   ?sub a schema:DataCatalog .
#'   ?subtype a schema:DataType .
#' }"
#' result_df <- .querySPARQL(endpoint, query)
#'
#' @importFrom httr content GET timeout add_headers
#' @importFrom utils URLencode read.csv
#' @noRd
.querySPARQL <- function(query) {

    result <- httr::content(
        httr::GET(
            query,
            httr::timeout(60),
            httr::add_headers(c(Accept = "text/csv"))
        ),
        "text", encoding = "UTF-8"
    )

    utils::read.csv(textConnection(result), stringsAsFactors = TRUE)
}



#
# prefixes <-
#
#     "PREFIX up_core: <http://purl.uniprot.org/core/>
# PREFIX uniref: <http://purl.uniprot.org/uniref/>
# PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
# "
#
# # May not be used
# data.frame(
#     name = c("species", "uniref90", "uniref50", "ec"),
#     code = c("sciName", "unirefId", "unirefId", "ec"),
#     pred = c("up_core:scientificName", "uniref:unirefId",
#              "uniref:unirefId", "up_core:enzyme")
# )
#
#
#
# # x <- levels(
# #     ariadne::importMapping(
# #         ariadne::ariadne("choco", uniref90 ~ ko),
# #         dry.run = FALSE
# #         )
# #     )[["uniref90"]][1:100]
#
#
# x <-
#     c("UniRef90_A0A010PZR5", "UniRef90_A0A010PZT4", "UniRef90_A0A010PZU0",
#       "UniRef90_A0A010PZV8", "UniRef90_A0A010PZW7", "UniRef90_A0A010Q006",
#       "UniRef90_A0A010Q047", "UniRef90_A0A010Q0B1", "UniRef90_A0A010Q0B7",
#       "UniRef90_A0A010Q0D0", "UniRef90_A0A010Q0G4", "UniRef90_A0A010Q0J0",
#       "UniRef90_A0A010Q0K3", "UniRef90_A0A010Q0U1", "UniRef90_A0A010Q136",
#       "UniRef90_A0A010Q165", "UniRef90_A0A010Q1D5", "UniRef90_A0A010Q1G0",
#       "UniRef90_A0A010Q1I2", "UniRef90_A0A010Q1I9", "UniRef90_A0A010Q1R3",
#       "UniRef90_A0A010Q1S8", "UniRef90_A0A010Q1W2", "UniRef90_A0A010Q216",
#       "UniRef90_A0A010Q273", "UniRef90_A0A010Q2D9", "UniRef90_A0A010Q2E2",
#       "UniRef90_A0A010Q2I5", "UniRef90_A0A010Q2J6", "UniRef90_A0A010Q2N4",
#       "UniRef90_A0A010Q2R1", "UniRef90_A0A010Q3J8", "UniRef90_A0A010Q3K6",
#       "UniRef90_A0A010Q3S3", "UniRef90_A0A010Q3U3", "UniRef90_A0A010Q3W2",
#       "UniRef90_A0A010Q3X9", "UniRef90_A0A010Q454", "UniRef90_A0A010Q484",
#       "UniRef90_A0A010Q4B9", "UniRef90_A0A010Q4D5", "UniRef90_A0A010Q4J6",
#       "UniRef90_A0A010Q4L3", "UniRef90_A0A010Q4P4", "UniRef90_A0A010Q4T7",
#       "UniRef90_A0A010Q519", "UniRef90_A0A010Q536", "UniRef90_A0A010Q554",
#       "UniRef90_A0A010Q578", "UniRef90_A0A010Q5E7", "UniRef90_A0A010Q5Q2",
#       "UniRef90_A0A010Q5R3", "UniRef90_A0A010Q5X4", "UniRef90_A0A010Q5X9",
#       "UniRef90_A0A010Q5Y3", "UniRef90_A0A010Q5Y8", "UniRef90_A0A010Q644",
#       "UniRef90_A0A010Q688", "UniRef90_A0A010Q6H4", "UniRef90_A0A010Q6L2",
#       "UniRef90_A0A010Q6L8", "UniRef90_A0A010Q6N7", "UniRef90_A0A010Q6R4",
#       "UniRef90_A0A010Q6S9", "UniRef90_A0A010Q732", "UniRef90_A0A010Q7B7",
#       "UniRef90_A0A010Q7D6", "UniRef90_A0A010Q7M8", "UniRef90_A0A010Q805",
#       "UniRef90_A0A010Q813", "UniRef90_A0A010Q832", "UniRef90_A0A010Q859",
#       "UniRef90_A0A010Q892", "UniRef90_A0A010Q8B0", "UniRef90_A0A010Q8B7",
#       "UniRef90_A0A010Q8E9", "UniRef90_A0A010Q8F2", "UniRef90_A0A010Q8F7",
#       "UniRef90_A0A010Q8G6", "UniRef90_A0A010Q8H5", "UniRef90_A0A010Q8H9",
#       "UniRef90_A0A010Q8J3", "UniRef90_A0A010Q8N5", "UniRef90_A0A010Q8P2",
#       "UniRef90_A0A010Q8P3", "UniRef90_A0A010Q8U1", "UniRef90_A0A010Q8Y7",
#       "UniRef90_A0A010Q8Y9", "UniRef90_A0A010Q8Z8", "UniRef90_A0A010Q910",
#       "UniRef90_A0A010Q922", "UniRef90_A0A010Q984", "UniRef90_A0A010Q985",
#       "UniRef90_A0A010Q9A4", "UniRef90_A0A010Q9B7", "UniRef90_A0A010Q9E1",
#       "UniRef90_A0A010Q9E5", "UniRef90_A0A010Q9E7", "UniRef90_A0A010Q9H0",
#       "UniRef90_A0A010Q9I7")
#
#
# sparkmap(x, endpoint = "uniprot", prefixes = prefixes)
#
