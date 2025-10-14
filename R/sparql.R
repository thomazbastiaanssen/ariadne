# The following functions are from https://github.com/aourednik/SPARQLchunks
# under the GPL-3 licence.
# We will temporarily mirror them here, eagerly awaiting their  ongoing CRAN
# submission.

#' Fetch data from a SPARQL endpoint and store the output in a dataframe
#' @param endpoint The SPARQL endpoint (a URL)
#' @param query The SPARQL query (character)
#' @param autoproxy Try to detect a proxy automatically (boolean). Useful on Windows machines behind corporate firewalls
#' @param auth Authentication Information (httr-authenticate-object)
#' @return SPARQL query result in data.frame format
#' @examples
#' endpoint <- "https://lindas.admin.ch/query"
#' query <- "PREFIX schema: <http://schema.org/>
#'   SELECT * WHERE {
#'   ?sub a schema:DataCatalog .
#'   ?subtype a schema:DataType .
#' }"
#' result_df <- sparql2df(endpoint, query)
#' @export
sparql2df <- function(endpoint, query, autoproxy = FALSE, auth = NULL) {
    proxy_config <- ifelse(
        autoproxy,
        autoproxyconfig(endpoint),
        httr::use_proxy(url = NULL)
    )
    acceptype <- "text/csv"
    outcontent <- get_outcontent(endpoint, query, acceptype, proxy_config, auth)
    tryCatch(
        content <- textConnection(outcontent$content),
        error = function(e) {
            stop(
                sprintf(
                    "There is something wrong with the output content: %s",
                    outcontent
                )
            )
        }
    )
    tryCatch(
        {
            df <- utils::read.csv(content)
        },
        error = function(e) {
            # utils::browseURL(outcontent$httpquery)
            stop(
                sprintf(
                    "Reply from SPARQL endpoint received but could not convert it to a data.frame.\nVerify the query result in a web browser:\n%s",
                    outcontent$httpquery
                )
            )
        }
    )
    return(df)
}

#' Get the content from the endpoint
#' @param endpoint The SPARQL endpoint (URL)
#' @param query The SPARQL query (character)
#' @param acceptype 'text/csv' or 'text/xml' (character)
#' @param proxy_config Detected proxy configuration (list)
#' @param auth Authentication Information (httr-authenticate-object)
#' @return The result of the SPARQL query as a list or, if this fails, failure message.
#' @noRd
get_outcontent <- function(endpoint, query, acceptype, proxy_config, auth = NULL) {
    qm <- paste(endpoint, "?", "query", "=",
                gsub("\\+", "%2B", utils::URLencode(query, reserved = TRUE)), "",
                sep = ""
    )
    message("Preparing to send query to: ", endpoint)
    message("SPARQL string:\n", query)
    message("Query URL:\n", qm )
    content <- tryCatch(
        {
            out <- httr::GET(
                qm,
                proxy_config, auth,
                httr::timeout(60),
                httr::add_headers(c(Accept = acceptype)),
                httr::user_agent("R client SPARQLChunks")
            )
            if (out$status == 401) {warning(
                "Authentication required. Provide valid authentication with the auth parameter"
            )} else {
                httr::warn_for_status(out)
            }
            httr::content(out, "text", encoding = "UTF-8") # Don't use return(...). If you use return(...) inside a block that is being assigned (x <- tryCatch({...})), you're exiting the function, not just returning a value for assignment.
        },
        error = function(e) {
            # @see https://github.com/r-lib/httr/issues/417
            # The download.file function in base R uses IE settings, including proxy password, when you use download
            # method wininet which is now the default on windows.
            if (is_windows()) {
                tempfile <- file.path(tempdir(), "temp.txt")
                utils::download.file(qm,
                                     method = "wininet",
                                     headers = c(Accept = acceptype),
                                     tempfile
                )
                temp <- paste(readLines(tempfile), collapse = "\n")
                unlink(tempfile)
                temp
            }
        }
    )
    if (is.null(content) || nchar(content) < 1) {
        warning(
            sprintf("First query attempt result is empty. Trying without Accept=%s header. The result is not guaranteed to be a list.",
                    acceptype)
        )
        content <- tryCatch(
            {
                out <- httr::GET(
                    qm,
                    proxy_config, auth,
                    httr::timeout(60),
                    httr::user_agent("R client SPARQLChunks")
                )
                if (out$status == 401) {warning(
                    "Authentication required. Provide valid authentication with the auth parameter"
                )} else {
                    httr::warn_for_status(out)
                }
                httr::content(out, "text", encoding = "UTF-8")
            },
            error = function(e) {
                if (is_windows()) {
                    tempfile <- file.path(tempdir(), "temp.txt")
                    utils::download.file(qm, method = "wininet", tempfile)
                    temp <- paste(readLines(tempfile), collapse = "\n")
                    unlink(tempfile)
                    temp
                }
            }
        )
        if (is.null(content) || nchar(content) < 1) {
            warning("The query result is still empty")
        }
    }
    if (inherits(content, "response")) {
        if (httr::status_code(content) >= 400) {
            stop(sprintf("HTTP error %s: %s", httr::status_code(content), httr::http_status(content)$message))
        }
    }
    return(list(
        content = content,
        httpquery = qm
    ))
}

#' Try to determine the proxy settings automatically
#' @param endpoint The SPARQL endpoint (URL)
#' @return Confirmation of the proxy setting
#' @noRd
autoproxyconfig <- function(endpoint) {
    message("Trying to determine proxy parameters")
    proxy_url <- tryCatch(
        {
            curl::ie_get_proxy_for_url(endpoint)
        },
        error = function(e) {
            message("Automatic proxy detection with curl::curl::ie_get_proxy_for_url() failed.")
            return(NULL)
        }
    )
    if (!is.null(proxy_url)) {
        message(paste("Using proxy:", proxy_url))
    } else {
        message(paste("No proxy found, nor needed, to access the endpoint", endpoint))
    }
    return(httr::use_proxy(url = proxy_url))
}

#' Verify if platform is Windows
#' @return TRUE if platform is Windows, FALSE otherwise
#' @noRd
is_windows <- function() {
    .Platform$OS.type == "windows"
}



# Query taxonomies based on uniref ids from UniProt using SPARQL
# x = Character vector of uniref IDs
.querySPARQL <- function(x, endpoint = "https://sparql.uniprot.org/") {
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("uniref:", x, collapse = " ")
    # Define first part of query
    query_part1 <- "
        PREFIX uniprot: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        # Outer query to get final names
        SELECT ?name
        WHERE {
        {
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
   sparql2df(endpoint, query)
}
