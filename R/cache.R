
#' @importFrom BiocFileCache bfcquery
#' @importFrom arrow write_parquet
.cache_resource <- function(url, res.name, from, to){
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
            x, key.FUN = function(keys) sub("GO:", "", keys, fixed = TRUE)
        ),
        WoL = function(x) .process_one2many(
            x, key.FUN = ifelse(from == "uniref90",
                function(keys) paste0("UniRef90_", keys), identity),
            val.FUN = function(vals) sub("EC-", "", vals, fixed = TRUE)
        ),
        BugSigDB = function(x) .process_one2many(
            x, val.cols = -c(1L, 2L), skip = 1L, key.FUN = function(keys){
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
        GM = .process_complex_modules,
        MSigDB = .process_rdslist
    )
    # Preprocess data
    linkmap <- FUN(url)
    # Add colnames
    if( res.name != "MSigDB" ) colnames(linkmap) <- c(from, to)
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
#' @importFrom stringr str_split fixed
.process_one2many <- function(x, key.col = 1L, val.cols = -key.col,
    key.FUN = identity, val.FUN = identity,...){
    # Read file content
    x <- read_lines(x, ...)
    # Split elements in each line by tab
    line_content <- str_split(x, fixed("\t"))
    # Extract keys
    keys <- vapply(line_content, `[`, key.col, FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line_content, `[`, val.cols)
    # Apply custom processing function to keys
    keys <- key.FUN(keys)
    # Create linkmap
    linkmap <- data.frame(
        x = rep(keys, lengths(values, use.names = FALSE)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    # Apply custom processing function to values
    linkmap$y <- val.FUN(linkmap$y)
    return(linkmap)
}


#' @importFrom BiocParallel bplapply
#' @importFrom data.table rbindlist
#' @importFrom utils download.file unzip
.process_rdslist <- function(x){
    # Prepare temporary file and directory
    temp_file <- tempfile(fileext = ".zip")
    temp_dir <- tempfile("dir")
    # Download file archive from Zenodo resource
    download.file(x, temp_file)
    unzip(temp_file, exdir = temp_dir)
    # List files from zipped file
    files <- list.files(temp_dir, full.names = TRUE)
    # Unzip file if compressed
    if( length(files) == 1L && endsWith(files, ".zip") ){
        unzip(files, exdir = temp_dir)
    }
    # List relevant unzipped files
    files <- list.files(temp_dir, full.names = TRUE)
    files <- grepv("zip|summary", files, invert = TRUE)
    # Import and bind linkmaps
    linkmaps <- bplapply(files, readRDS)
    linkmap <- rbindlist(linkmaps)
    # Remove temporary file and directory
    unlink(c(temp_file, temp_dir), recursive = TRUE)
    return(linkmap)
}
