
#' @export
importMapping <- function(map.file){
  
    # Read the file content
    lines <- readLines(map.file)
  
    values <- list()
    keys <- list()
  
    for (line in lines) {
        # Split by tab character
        items <- unlist(strsplit(line, "\t"))
        # Append the entire line as a vector
        values <- append(vectors, list(items[-1]))
        keys <- append(keys, items[1])
  }
  
  names(values) <- keys
  return(values)
}

module.bases <- list(
    GMM = "https://raw.githubusercontent.com/raeslab/GMMs/refs/heads/master/GMMs.v1.07.txt",
    # GBM = "https://raeslab.org/software/GBMs.zip"
)

#' @export
importModules <- function(module.file){
  
    if( module.file %in% names(module.bases) ){
        lines <- readLines(module.bases[[module.file]])
    }else{
        lines <- readLines(module.file)
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