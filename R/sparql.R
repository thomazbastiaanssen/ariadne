#' @title SPARQLmap
#' @name SPARQLmap
#' @param x `Character vector`. Optional. Feature IDs on which the query should
#'     be constrained.
#' @param x.type `Character scalar`. Feature type of `x`. (Default: `"uniref"`)
#' @param .by either a `formula`, or a `Character vector` of length 2,  with the
#'      names of the desired combination of feature types.
#' @param endpoint The SPARQL endpoint (a URL)
#' @returns a data.frame of desired query.
#' @export
#' @examples
#'
#'
#'x <- paste0(
#'    "UniRef90_",
#'    c("A0A010PZR5", "A0A010PZT4", "A0A010PZU0", "A0A010PZV8",
#'    "A0A010PZW7", "A0A010Q006", "A0A010Q047", "A0A010Q0B1")
#'    )
#'
#'
#' SPARQLmap(x, .by = uniref ~ ec,  endpoint = "uniprot")
#'
#' SPARQLmap(x, .by = uniref ~ species,  endpoint = "uniprot")
#'
SPARQLmap <- function(
        x, x.type = "uniref", .by = NULL, endpoint = "uniprot"
        ) {

    terms <- if (inherits(.by, "formula"))
        all.vars(.by)
    else .by

    known_endpoints <- c("uniprot")
    endpoint <- known_endpoints[charmatch(endpoint, table = known_endpoints)]
    endpoint_url <- switch(endpoint, "uniprot" = "https://sparql.uniprot.org/")

    query <- .composeSPARQL(x, x.type, terms, endpoint_url)
    output <- .querySPARQL(query)
    colnames(output) <- terms
    return(output)
}

#' @importFrom utils URLencode
.composeSPARQL <- function(x, x.type, terms, endpoint_url) {
    #uniref.ids <- paste0("uniref:", x, collapse = " ")
 terms_remote <- vapply(terms, .uniprot_terms, "")

 prefixes <-
"PREFIX up_core: <http://purl.uniprot.org/core/>
PREFIX uniref: <http://purl.uniprot.org/uniref/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
"

 uniref2taxonomy <-
"?unirefId up_core:member ?member .
?member up_core:organism ?taxId .
?taxId up_core:scientificName ?sciName .
"

 uniref2ec <-
"?protein up_core:representativeFor ?unirefId .
?protein ( up_core:enzyme | up_core:domain/up_core:enzyme | up_core:component/up_core:enzyme ) ?ecId .
"

q <- paste0("
SELECT", paste0(terms_remote, collapse = " "), "
WHERE {
", .composeX(x, x.type), "
",
paste0(
    if("?ecId" %in% terms_remote) uniref2ec,
    if("?sciName" %in% terms_remote) uniref2taxonomy
    ),
"}")
    # Build query
    query <- paste0(prefixes, q)
    query <- gsub("\\+", "%2B", utils::URLencode(query, reserved = TRUE))
    query <- paste0(endpoint_url, "?query=", query)
    return(query)
}

.composeX <- function(x, x.type) {
    if(length(x) == 0L) return("")
    if(x.type %in% c("uniref", "uniref50", "uniref90", "uniref100")) {
        return(
            paste0(
                "VALUES ?unirefId {
            ", paste0("uniref:", x, collapse = " "), "
        } .")
        )
    }
    if(x.type %in% c("ec")) {
        return(
            paste0(
                "VALUES ?ecId {
            ", paste0("enzyme:", x, collapse = " "), "
        } .")
        )
    }
    if(x.type %in% c("species")) {
        return(
            paste0(
                "VALUES ?SciName {
            ", paste0("up_core:", x, collapse = " "), "
        } .")
        )
    }
}

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

.uniprot_terms <- function(x) {
    switch(x,
           "uniref"    = ,
           "uniref50"  =,
           "uniref90"  =,
           "uniref100" = "?unirefId",
           "species"   = "?sciName",
           "ec"        = "?ecId")
}





