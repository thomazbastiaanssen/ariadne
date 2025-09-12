
#' @export
importMapping <- function(map.file){
    # Whether to use package or custom mapping
    if( map.file %in% names(ChocoPhlAn) ){
        # Derive col.names
        col.names <- unlist(strsplit(gsub(
            ".*_([^_]+_[^_]+)\\.txt\\.gz$", "\\1", ChocoPhlAn[[map.file]]), "_")
        )
        # Cache database
        map.file <- .getFile(ChocoPhlAn[[map.file]])
    }
    # Read file content
    lines <- readLines(map.file)
    # Split elements in each line by tab
    items <- strsplit(x = lines, split = "\t")
    # Extract keys
    keys <- vapply(items, FUN = function(x) x[1], FUN.VALUE = "")
    # Extract values
    values <- lapply(items, FUN = function(x) x[-1])
    # Import as linkmap
    linkmap <- as.linkMap(keys, values, col.names)
    return(linkmap)
}

#' @export
importModules <- function(module.file){
    # Whether to use package or custom modules
    if( module.file %in% names(GutModules) ){
        # Cache database
        module.file <- .getFile(GutModules[[module.file]])
    }
    # Read the file content
    lines <- readLines(module.file)
    # GBM is missing /// at the end
    if( module.file == "GBM" ){
        lines <- append(lines, "///")
    }
    
    keys <- c()
    modules <- list()
    module <- c()
    
    for (line in lines){
      
        if( line == "///" ){
        
            module <- unlist(strsplit(module, ",|\t"))
            modules <- append(modules, list(module))
            module <- c()
        
        }else if( substr(line, 1, 1) == "M" ){
        
            key <- paste0(gsub("\t", " (", line), ")")
            keys <- append(keys, key)
        
        }else{
        
            module <- append(module, line)
        
        }
    }
    
    linkmap <- as.linkMap(keys, modules, col.names = c("mod", "ko"))
    return(linkmap)
}

#' @export
as.linkMap <- function(keys, values, col.names = NULL){
    # Create linkMap
    linkmap <- data.frame(
        x = rep(keys, vapply(values, length, 1)),
        y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
    # Assign custom colnames
    if( !is.null(col.names) ){
        names(linkmap) <- col.names
    }
    return(linkmap)
}
