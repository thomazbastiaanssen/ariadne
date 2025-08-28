

#' @export
importMapping <- function(map.file){
  
    if( map.file %in% names(mapDatabases) ){
        # Cache database
        cached <- .getFile(mapDatabases[[map.file]])
        # Read the file content
        lines <- readLines(cached)
    }else{
        # Read the file content
        lines <- readLines(map.file)
    }
    
    values <- list()
    keys <- c()
  
    for (line in lines) {
        # Split by tab character
        items <- unlist(strsplit(line, "\t"))
        # Append the entire line as a vector
        values <- append(values, list(items[-1]))
        keys <- append(keys, items[1])
    }
  
    names(values) <- keys
  
    values <- lapply(values, function(x) gsub("UniRef90_", "", x))
  
    return(values)
}

#' @export
importModules <- function(module.file){
  
    if( module.file %in% names(moduleDatabases) ){
    
        cached <- .getFile(moduleDatabases[[module.file]])
        lines <- readLines(cached)
    
    }else{
    
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
