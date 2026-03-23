#' Utility functions
#' @name utils
#' @rdname utils
#' 
#' @description
#' These utility functions are used throughout the package and may be relevant
#' in other packages dealing with annotation mappings. \code{as.linkmap}
#' converts a list of named vectors to a linkmap data.frame.
#' 
#' @param se A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @param x The rowData or colData of a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @param modules \code{data.frame}. A linkmap returned by weavePath or
#'   weaveComplex.
#' 
#' @param by \code{Character scalar} A string specifying a variable of \code{x}
#'   to append \code{modules}. (Default: \code{"row.names"})
#' 
#' @returns
#' \code{as.linkmap} returns a linkmap \code{data.frame} where the first and
#' second columns contains \code{keys} and \code{values} and each row represents
#' a unique combination of the two.
NULL


#' @export
#' @rdname utils
#' @importFrom SummarizedExperiment rowData
#' @importFrom stringr fixed str_detect str_split
processGeneFamilies <- function(se){
    # Select rows with non-null taxa
    se <- se[str_detect(rownames(se), fixed("|")), ]
    se <- se[str_detect(rownames(se), "unclassified", negate = TRUE), ]
    # Split gene and taxonomy
    gene.linkmap <- as.data.frame(
        str_split(rownames(se), fixed("|"), n = 2, simplify = TRUE)
    )
    names(gene.linkmap) <- c("uniref90", "taxname")
    # Split genus and species
    tax.linkmap <- as.data.frame(
        str_split(gene.linkmap$taxname, fixed("."), n = 2, simplify = TRUE),
    )
    names(tax.linkmap) <- c("genus", "species")
    # Bind gene and tax linkmaps
    rowData(se) <- cbind(rowData(se), gene.linkmap, tax.linkmap)
    return(se)
}


#' @export
#' @rdname utils
appendModules <- function(x, modules, by = "row.names"){
    # Check args
    if( !by %in% c("row.names", colnames(x)) ){
        stop("'by' must be 'row.names' or a variable of 'x'.", call. = FALSE)
    }
    # Convert linkmap to wide format
    modules <- table(modules) == 1
    
    if( by == "row.names" ){
        idx <- match(rownames(x), rownames(modules))
    }else{
        idx <- match(x[[by]], rownames(modules))
    }
    
    modules <- modules[idx, , drop = FALSE]
    modules[is.na(modules)] <- FALSE
    
    out <- cbind(x, modules)
    return(out)
}
