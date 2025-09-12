#' Utility functions
#'
#' @name utils
NULL

#' @rdname utils
#' @export
setMethod("as.linkMap", signature = c(values = "list"),
    function(values, keys = NULL, col.names = NULL){
        if( is.null(keys) ){
            keys <- names(values)
        }
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
)