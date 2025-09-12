#' Import modules from a file or a database
#' 
#' @name importModules
NULL

ModuleDatabases <- list(
    GMM = "GMMs.v1.07.txt",
    GBM = "GBMs.v1.0.txt"
)

#' @rdname importModules
#' @export
setMethod("importModules", signature = c(module.file = "character"),
    function(module.file){
        # Whether to use package or custom modules
        if( module.file %in% names(ModuleDatabases) ){
            # Cache database
            module.file <- .getCache(ModuleDatabases[[module.file]])
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
)