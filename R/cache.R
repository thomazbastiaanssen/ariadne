

#' @importFrom BiocFileCache bfcquery
#' @importFrom utils download.file
#' @importFrom arrow write_parquet
.cache_resource <- function(url, res.name, from, to) {
    # Initialise cache
    bfc <- .init_cache()
    # Build resource name
    rname <- file.path(res.name, basename(url))
    # Check if preprocessed file is cached
    cached <- bfcquery(bfc, rname)$rpath
    # Return preprocessed file if available
    if( length(cached) > 0L ) return(cached[1])
    # Select function based on resource
    FUN <- switch(
        res.name,
        ChocoPhlAn = function(x) .process_one2many(
            x, FUN = function(keys) sub("GO:", "", keys, fixed = TRUE)
        ),
        WoL = function(x) {
            # Download temporary file (xz not supported by read_lines)
            temp_xz <- tempfile(fileext = ".xz")
            download.file(x, temp_xz, mode = "wb", quiet = TRUE)
            # Process temporary file
            .process_one2many(
                temp_xz, FUN = function(keys) paste0("UniRef90_", keys)
            )
        },
        BugSigDB = function(x) .process_one2many(
            x, nonval.cols = c(1L, 2L), skip = 1L, FUN = function(keys){
                # Remove module prefix
                keys <- sub("bsdb:", "", keys, fixed = TRUE)
                # Remove module description
                keys <- sub("_.*$", "", keys)
            }
        ),
        TIGRFAMs = function(x) .process_one2one(
            x, header = FALSE, select = c(1L, 2L)
        ),
        GO = function(x) .process_one2one(x, header = FALSE),
        Misc = .process_complex_modules
    )
    # Preprocess data
    linkmap <- FUN(url)
    # Add colnames
    colnames(linkmap) <- c(from, to)
    # Store linkmap in cache as parquet file
    .add2cache(linkmap, rname, bfc)
    # Return path to cached preprocessed file
    cached <- bfcquery(bfc, rname)$rpath
    return(cached)
}


#' @importFrom BiocFileCache BiocFileCache
#' @importFrom tools R_user_dir
.init_cache <- function(){
    # Define cache dir
    cache <- R_user_dir("ariadne", "cache")
    # Initialise cache
    bfc <- BiocFileCache(cache, ask = FALSE)
    return(bfc)
}


#' @importFrom arrow write_parquet
#' @importFrom BiocFileCache bfccache bfcadd
.add2cache <- function(x, rname, bfc){
    # Find path for cache subdir
    fpath <- bfc |>
        bfccache() |>
        file.path(rname)
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
    linkmap <- fread(x, ...)
    # Remove id prefix ending with : (for GO resources)
    linkmap$V1 <- sub("^[^:]*:", "", linkmap$V1)
    # Remove GO prefix (for GO and TIGRFAMs resources)
    linkmap$V2 <- sub("GO:", "", linkmap$V2, fixed = TRUE)
    return(linkmap)
}


#' @importFrom readr read_lines
.process_one2many <- function(
    x, key.col = 1L, nonval.cols = key.col, FUN = identity, ...){
    # Read file content
    x <- read_lines(x, ...)
    # Split elements in each line by tab
    line_content <- strsplit(x, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line_content, `[`, key.col, FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line_content, `[`, -nonval.cols)
    # Apply custom processing function
    keys <- FUN(keys)
    # Create linkmap
    linkmap <- data.frame(
        x = rep(keys, lengths(values, use.names = FALSE)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    return(linkmap)
}
