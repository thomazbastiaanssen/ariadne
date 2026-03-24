

#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd
#' @importFrom tools R_user_dir
#' @importFrom arrow write_parquet
.cache_resource <- function(url, resource, from, to) {
    # Initialise cache
    cache <- R_user_dir("ariadne", "cache")
    bfc <- BiocFileCache(cache, ask = FALSE)
    # Build resource name
    rname <- file.path(resource, basename(url))
    # Check if preprocessed file is cached
    cached <- bfcquery(bfc, rname)$rpath
    # Return preprocessed file if available
    if (length(cached) > 0L) return(cached[1])
    # Find path for cache subdir
    path <- file.path(cache, rname)
    subdir <- dirname(path)
    # Create subdir if absent
    if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
    # Download raw file
    download.file(url, path)
    # Select function based on resource
    FUN <- switch(
        resource,
        BugSigDB = .process_bugsigdb,
        ChocoPhlAn = .process_chocophlan,
        GM = .process_complex_modules,
        GO = .process_go,
        TIGRFAMs = .process_tigrfams,
        WoL = .process_wol
    )
    # Read file content
    x <- readLines(path)
    # Preprocess data
    linkmap <- FUN(x)
    # Add colnames
    colnames(linkmap) <- c(from, to)
    # Store preprocessed data
    write_parquet(linkmap, path)
    # Cache the preprocessed file using URL as resource name
    bfcadd(
        bfc, rname = rname, fpath = path, fname = "exact", action = "asis"
    )
    # Return path to cached preprocessed file
    cached <- bfcquery(bfc, rname)$rpath
    return(cached)
}


.process_chocophlan <- function(x){
    # Split elements in each line by tab
    line.content <- strsplit(x, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line.content, `[`, -1L)
    # Create linkmap
    linkmap <- data.frame(
        x = rep(keys, lengths(values, use.names = FALSE)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    # Remove GO id prefix ending with :
    linkmap$x <- gsub("GO:", "", linkmap$x, fixed = TRUE)
    return(linkmap)
}


.process_wol <- function(x){
    linkmap <- .process_chocophlan(x)
    linkmap$x <- paste0("UniRef90_", linkmap$x)
    return(linkmap)
}


.process_go <- function(x){
    # Remove header
    x <- x[!startsWith(x, "!")]
    # Split entries into keys and values
    x <- sub("^(\\S+).*?(\\S+)$", "\\1 \\2", x)
    x <- strsplit(x, " ", fixed = TRUE)
    # Create linkmap
    linkmap <- as.data.frame(do.call(rbind, x))
    # Trim prefix ending with :
    linkmap <- data.frame(
        x = gsub("^[^:]*:", "", linkmap[, 1L]),
        y = gsub("^[^:]*:", "", linkmap[, 2L])
    )
    return(linkmap)
}


.process_tigrfams <- function(x){
    # Split elements in each line by tab
    line.content <- strsplit(x, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = character(1L))
    # Extract values
    values <- vapply(line.content, `[`, 2L, FUN.VALUE = character(1L))
    # Remove id prefix ending with :
    values <- gsub("^[^:]*:", "", values)
    # Create linkmap
    linkmap <- data.frame(x = keys, y = values)
    return(linkmap)
}


.process_bugsigdb <- function(x){
    # Remove header
    x <- x[-1L]
    # Split elements in each line by tab
    line.content <- strsplit(x, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = character(1L))
    # Remove module prefix
    keys <- sub("bsdb:", "", keys, fixed = TRUE)
    # Remove module description
    keys <- sub("_.*$", "", keys)
    # Extract values
    values <- lapply(line.content, `[`, -c(1L, 2L))
    # Create linkmap
    linkmap <- data.frame(
        x = rep(keys, lengths(values, use.names = FALSE)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    return(linkmap)
}
