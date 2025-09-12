
#' @export
#' @importFrom utils read.table
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
    # Add names to list
    names(values) <- keys
    return(values)
}

importMappings <- function(map.files, merge = TRUE){
    maps <- lapply(map.files, importMapping)
    if( merge ){
        maps <- unlist(maps, recursive = FALSE)
    }
    return(maps)
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
    if( lines[length(lines)] != "///" ){
        lines <- append(lines, "///")
    }
    
    keys <- c()
    modules <- list()
    module <- c()
    
    for (line in lines){
      
        if( line == "///" ){
        
            modules <- append(modules, list(module))
            module <- c()
        
        }else if( substr(line, 1, 1) == "M" ){
        
            key <- paste0(gsub("\t", " (", line), ")")
            keys <- append(keys, key)
        
        }else{
            module <- append(module, line)
        }
    }
    
    for( i in seq_along(modules) ){
        module.comp <- strsplit(modules[[i]], "\t")
        modules[[i]] <- lapply(module.comp, strsplit, split = ",")
    }
    
    names(modules) <- keys
    return(modules)
}
