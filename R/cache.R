

#' @importFrom BiocFileCache BiocFileCache bfcquery
#' @importFrom tools R_user_dir
#' @importFrom utils download.file
#' @importFrom arrow write_parquet
.cache_resource <- function(url, res.name, from, to) {
    # Initialise cache
    cache <- R_user_dir("ariadne", "cache")
    bfc <- BiocFileCache(cache, ask = FALSE)
    # Build resource name
    rname <- file.path(res.name, basename(url))
    # Check if preprocessed file is cached
    cached <- bfcquery(bfc, rname)$rpath
    # Return preprocessed file if available
    if( length(cached) > 0L ) return(cached[1])
    # Find path for cache subdir
    path <- file.path(cache, rname)
    # Download raw file
    download.file(url, path)
    # Select function based on resource
    FUN <- switch(
        res.name,
        BugSigDB = .process_bugsigdb,
        ChocoPhlAn = .process_chocophlan,
        GM = .process_complex_modules,
        WoL = .process_wol,
        TIGRFAMs = ,
        GO = function(x) process_one2one(x, header = FALSE)
    )
    # Read file content
    x <- readLines(path)
    # Preprocess data
    linkmap <- FUN(x)
    # Add colnames
    colnames(linkmap) <- c(from, to)
    # Store linkmap in cache as parquet file
    .add2cache(linkmap, rname, path, bfc)
    # Return path to cached preprocessed file
    cached <- bfcquery(bfc, rname)$rpath
    return(cached)
}


#' @importFrom arrow write_parquet
#' @importFrom BiocFileCache bfcadd
.add2cache <- function(x, rname, fpath, bfc){
    # Define resource subdir of cache dir
    subdir <- dirname(fpath)
    # Create subdir if absent
    if( !dir.exists(subdir) ) dir.create(subdir, recursive = TRUE)
    # Store preprocessed data
    write_parquet(x, fpath)
    # Cache the preprocessed file using URL as resource name
    bfcadd(
        bfc, rname = rname, fpath = fpath, fname = "exact", action = "asis"
    )
}


# from to args?
#' @importFrom data.table fread
.process_one2one <- function(x, ...){
    # Read linkmap
    linkmap <- fread(text = x, ...)
    # Remove id prefix ending with : (for GO resources)
    linkmap$V1 <- sub("^[^:]*:", "", linkmap$V1)
    # Remove GO prefix (for GO and TIGRFAMs resources)
    linkmap$V2 <- sub("GO:", "", linkmap$V2, fixed = TRUE)
    # Select appropriate columns (for TIGRFAMs resources)
    linkmap <- linkmap[ , c(1L, 2L)]
    return(linkmap)
}


.process_chocophlan <- function(x){
    # Split elements in each line by tab
    line.content <- strsplit(x, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line.content, `[`, -1L)
    # Remove GO id prefix ending with :
    keys <- gsub("GO:", "", keys, fixed = TRUE)
    # Create linkmap
    linkmap <- data.frame(
        x = rep(keys, lengths(values, use.names = FALSE)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    return(linkmap)
}


.process_wol <- function(x){
    linkmap <- .process_chocophlan(x)
    linkmap$x <- paste0("UniRef90_", linkmap$x)
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
