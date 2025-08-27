

# Adapted from UniProt.ws utilities.R
#' @importFrom BiocFileCache BiocFileCache bfcneedsupdate bfcrpath bfcdownload
.getFile <- function(url){
  
    cache <- tools::R_user_dir("ariadne", "cache")
    
    bfc <- BiocFileCache(cache, ask = FALSE)
    
    rpath <- bfcrpath(
        bfc, rnames = url, exact = TRUE, download = TRUE, rtype = "web"
    )
    
    update <- bfcneedsupdate(bfc, names(rpath))
    
    if (update){
        bfcdownload(bfc, names(rpath), ask = FALSE)
    }

    return(rpath)
}
