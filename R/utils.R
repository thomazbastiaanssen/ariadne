#' Utility functions
#' 
#' These utility functions are used throughout the package and may be relevant
#' in other packages dealing with annotation mappings. \code{as.linkmap}
#' converts a list of named vectors to a linkmap data.frame.
#' 
#' @param values \code{Character list}. List of vectors where each vector
#'   represents one type of values.
#' 
#' @param keys \code{Character vector}. The vector of keys to use as names for
#'   \code{values} if the latter is unnamed. (Default: \code{NULL}).
#' 
#' @param col.names \code{Character vector}. A vector of two elements specifying
#'   the names of the columns in the output linkmap (Default: \code{NULL}).
#' 
#' @returns
#' \code{as.linkmap} returns a linkmap \code{data.frame} where the first and
#' second columns contains \code{keys} and \code{values} and each row represents
#' a unique combination of the two.
#'
#' @name utils
NULL

#' @rdname utils
#' @export
setMethod("as.linkmap", signature = c(values = "list"),
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