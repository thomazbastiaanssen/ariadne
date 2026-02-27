

#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd
#' @importFrom tools R_user_dir
#' @importFrom arrow write_parquet
.cache_resource <- function(url, resource) {
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
        ChocoPhlAn = .process_chocophlan,
        GBM = .process_complex_modules,
        GMM = .process_complex_modules,
        GO = .process_go,
        TIGRFAMs = .process_tigrfams,
        WoL = .process_wol
    )
  
    # Read file content
    x <- readLines(path)
    # Preprocess data
    linkmap <- FUN(x)
    # Store preprocessed data
    write_parquet(linkmap, path)
    # Cache the preprocessed file using URL as resource name
    bfcadd(
        bfc,
        rname = rname,
        fpath = path,
        fname = "exact",
        action = "asis"
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
    x <- .process_chocophlan(x)
    x <- x[ , c(2, 1)]
    return(x)
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
        x = gsub("^[^:]*:", "", linkmap[, 1]),
        y = gsub("^[^:]*:", "", linkmap[, 2])
    )
    return(linkmap)
}

.process_complex_modules <- function(x, br = "///", AND = ",", OR = "\t") {
    # Identify break lines
    v_br <- x == br
    # Split content by breaks, excluding break lines themselves
    line.content <- split(x[!v_br], cumsum(v_br)[!v_br])
    # Extract keys (first line of each block)
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = "", USE.NAMES = FALSE)
    # Extract values (all lines except first in each block)
    values <- lapply(line.content, `[`, -1L)
    # Replace tabs with spaces in keys, repeated for each value line
    module <- gsub("\t.*", "", rep(keys, lengths(values, use.names = FALSE)))
    # Create unique module_component identifiers
    module_component <- paste0(
        module, "_part_",
        unlist(lapply(rle(module)$lengths, seq), use.names = FALSE)
    )
    # Split feature list by tab character
    feature_list <- strsplit(
        unlist(values, recursive = TRUE, use.names = FALSE), "\t"
    )

    names(feature_list) <- module_component
    # Flatten feature complex list and split by comma to get individual features
    feature_complex <- unlist(feature_list, use.names = FALSE)
    feature <- strsplit(feature_complex, ",")
    # Create and return structured list of data frames
    modules <- list(
        module2module_component = data.frame(module,module_component),
        component2feature_complex = data.frame(
            module_component = rep(module_component, lengths(feature_list)),
            feature_complex
        ),
        feature_complex2feature = data.frame(
            feature_complex = rep(feature_complex, lengths(feature)),
            feature = unlist(feature, use.names = FALSE)
        )
    )
    modules <- MultiFactor(modules)
    linkmap <- weave(modules, module ~ feature)
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
