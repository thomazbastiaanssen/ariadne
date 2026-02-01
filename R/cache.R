
library(BiocFileCache)

#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd 
.getCache <- function(url, resource) {

  cache <- tools::R_user_dir("ariadne", "cache")
  bfc <- BiocFileCache(cache, ask = FALSE)
  
  rname <- file.path(resource, basename(url))
  
  # Check if preprocessed file is cached
  cached <- bfcquery(bfc, rname)$rpath
  
  # If preprocessed file found in cache
  if( length(cached) > 0L ){
      # Return path
      return(cached[1])
  }
  
  path <- file.path(cache, rname)
  subdir <- dirname(path)
  
  if( !dir.exists(subdir) ){
      dir.create(subdir)
  }
  
  # Download raw file
  download.file(url, path)

  FUN <- switch(
      resource,
      ChocoPhlAn = .process_chocophlan,
      WoL = .process_wol,
      GO = .process_go,
      #GM = .process_complex_modules,
      # KEGG = maybe no need to cache?
  )

  # Read file content
  x <- readLines(path)
  # Preprocess data
  linkmap <- FUN(x)
  # Store preprocessed data
  write.csv(linkmap, path, row.names = FALSE)
  
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
    
    linkmap <- data.frame(
      x = rep(keys, lengths(values, use.names = FALSE)),
      y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    
    return(linkmap)
}

.process_wol <- function(x){
    x <- .process_chocophlan(x)
    x <- x[ , c(2, 1)]
    return(x)
}

.process_go <- function(x){
    x <- x[!startsWith(x, "!")]
    x <- sub("^(\\S+).*?(\\S+)$", "\\1 \\2", x)
    x <- strsplit(x, " ", fixed = TRUE)
    x <- as.data.frame(do.call(rbind, x))
    return(x)
}



