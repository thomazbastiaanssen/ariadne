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
            x = rep(keys, lengths(values)),
            y = unlist(values, recursive = TRUE, use.names = FALSE)
        )
        # Assign custom colnames
        if( !is.null(col.names) ){
            names(linkmap) <- col.names
        }
        return(linkmap)
    }
)

#' @export
#' @rdname utils
#' @importFrom stringr fixed str_detect str_split
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<- assays
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
processGeneFamilies <- function(se){
    # Select rows with non-null taxa
    se <- se[str_detect(rownames(se), fixed("|")), ]
    se <- se[str_detect(rownames(se), "unclassified", negate = TRUE), ]
    # Split gene and taxonomy
    gene.linkmap <- as.data.frame(
        str_split(rownames(se), fixed("|"), n = 2, simplify = TRUE)
    )
    names(gene.linkmap) <- c("GeneID", "Taxon")
    # Split genus and species
    tax.linkmap <- as.data.frame(
        str_split(gene.linkmap$Taxon, fixed("."), n = 2, simplify = TRUE),
    )
    names(tax.linkmap) <- c("Genus", "Species")
    # Convert SE to TreeSE
    tse <- TreeSummarizedExperiment(
        assays = assays(se),
        rowData = cbind(rowData(se), gene.linkmap, tax.linkmap),
        colData = colData(se)
    )
    return(tse)
}
