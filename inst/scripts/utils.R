
fetch_zenodo_versions <- function(record_id){
    # Build version url
    url <- paste0("https://zenodo.org/api/records/", record_id, "/versions")
    # Request the versions list
    resp <- request(url) |>
        req_perform() |>
        resp_body_json()
    # Extract versions info
    hits <- resp$hits$hits
    # Obtain version and doi pairs
    versions <- as.character(vapply(hits, `[[`, "id", FUN.VALUE = integer(1L)))
    names(versions) <- vapply(hits, function(x) x$metadata$version, character(1L))
    return(versions)
}