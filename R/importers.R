
#' @export
#' @importFrom utils read.table
importMapping <- function(map.file){
    # Whether to use package or custom mapping
    if( map.file %in% names(mapDatabases) ){
        # Cache database
        map.file <- .getFile(mapDatabases[[map.file]])
    }
    # Read lines
    map <- read.table(
        map.file,
        sep = "\t",
        row.names = 1,
        fill = TRUE,
        flush = TRUE
    )
    # Convert data.frame to list of named vectors
    map <- apply(map, 1L, function(key)
        key[1:match("", key, nomatch = ncol(map) + 1) - 1])
    return(map)
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
