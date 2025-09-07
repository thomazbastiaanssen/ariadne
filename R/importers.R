

#' @export
importMapping <- function(map.file){
    # Whether to use package or custom mapping
    if( map.file %in% names(mapDatabases) ){
        # Cache database
        cached <- .getFile(mapDatabases[[map.file]])
        # Read the file content
        lines <- readLines(cached)
    }else{
        # Read the file content
        lines <- readLines(map.file)
    }
    # Convert each line to vector with key in 1 and values in 2+ positions
    items <- bplapply(lines, function(line) {
        item <- unlist(strsplit(line, "\t"))
        item <- list(key = item[1], value = item[-1])
        return(item)
    })
    # Retrieve keys and values
    keys <- vapply(items, function(x) x$key, character(1))
    values <- lapply(items, function(x) x$value)
    # Name values with keys
    names(values) <- keys
    return(values)
}

#' @export
importModules <- function(module.file){
    # Whether to use package or custom modules
    if( module.file %in% names(moduleDatabases) ){
        # Cache database
        cached <- .getFile(moduleDatabases[[module.file]])
        # Read the file content
        lines <- readLines(cached)
    }else{
        # Read the file content
        lines <- readLines(module.file)
    }
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
    
    names(modules) <- keys
    return(modules)
}
