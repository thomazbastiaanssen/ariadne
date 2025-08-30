#' @export
#' @importFrom dplyr bind_rows
mapModules <- function(modules, map){
    
    keep <- names(map) %in% unique(unlist(modules))
    
    map <- map[keep]
    
    tax <- lapply(map, .querySPARQL)
    
    tax.df <- bind_rows(tax, .id = "source")
    
    sig.list <- list()
    
    for( i in seq_along(modules) ){
      
        module <- names(modules[i])
        keep <- tax.df$source %in% modules[[i]]
    
        if( any(keep) ){
            sig.list[[module]] <- tax.df$name
        }
    }
    
    return(sig.list)
}
