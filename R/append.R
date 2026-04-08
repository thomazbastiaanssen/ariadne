#' Append linkmaps to SummarizedExperiment side information
#' 
#' @name addModules
#' 
#' @description
#' getModules and addModules allow to retrieve or append a linkmap to the side
#' information of a SummarizedExperiment (SE) object. This makes ariadne
#' interoperable with SE-based data analysis.
#' 
#' @param x The rowData or colData of a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @param modules \code{data.frame}. A linkmap as returned by
#'   \code{\link{weavePath}} or \code{\link{weaveComplex}}. Its first and second
#'   columns must contain elements to match to \code{key} and the target
#'   modules, respectively. 
#' 
#' @param by \code{Character scalar} A string indicating whether to append
#'   \code{modules} to the \code{"rows"} or \code{"cols"} side information of
#'   \code{x}. (Default: \code{"rows"})
#' 
#' @param key \code{Character vector} A vector specifying one or more variables
#'   of \code{x} side information based on which \code{modules} should be
#'   appended. (Default: \code{"row.names"})
#' 
#' @param as \code{Character scalar} A string specifying whether \code{modules}
#'   ids or names should be used. For the latter, a third column with names
#'   must exist in \code{modules}. (Default: \code{"ids"})
#' 
#' @param ... Unused.
#' 
#' @returns
#' An object of the same type as \code{x} with additional columns in its side
#' information, each containing information on membership to a certain module.
#' 
#' @examples
#' library(mia)
#' library(miaViz)
#' 
#' # Import datasets
#' data("Tengeler2020", package = "mia")
#' data("butyrate", package = "ariadne")
#' 
#' # Rename experiment object
#' tse <- Tengeler2020
#' 
#' # Get butyrate-producer module membership
#' modules <- getModules(tse, butyrate, key = "Genus")
#' 
#' # Get modules based on multiple variables given in order of priority
#' modules <- getModules(tse, butyrate, key = c("Genus", "Family"))
#' 
#' # Add modules to experiment as names instead of ids
#' tse <- addModules(tse, butyrate, key = "Genus", as = "names")
#' 
#' # Generate relative abundance table
#' tse <- transformAssay(tse, method = "relabundance")
#' 
#' # Agglomerate features by membership to butyrate-producer module
#' mod.se <- agglomerateByModule(tse, by = "rows", group = "butyrate")
#' 
#' # Plot relative abundance of butyrate producers
#' plotAbundance(mod.se, assay.type = "relabundance")
NULL


#' @export
#' @rdname addModules
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
setMethod("addModules", signature = c(x = "SummarizedExperiment"),
    function(x, modules, by = "rows", key = "row.names", as = "ids"){
    # Check margin
    if( length(by) == 0L || !by %in% c("rows", "cols") ){
        stop("'by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    # Select side information based on margin
    df <- if( by == "rows" ) rowData(x) else colData(x)
    # Build modules table
    modules <- .get_modules(df, modules, key, as)
    # Find any duplicate columns
    to_remove <- which(colnames(df) %in% colnames(modules))
    # Replace duplicate columns
    if( length(to_remove) != 0L ){
        warning("Some columns were replaced.", call. = FALSE)
        df[to_remove] <- NULL
    }
    # Append modules to side information
    out <- cbind(df, modules)
    if( by == "rows" ) rowData(x) <- out else colData(x) <- out
    return(x)
})


#' @export
#' @rdname addModules
#' @importFrom SummarizedExperiment rowData colData
setMethod("getModules", signature = c(x = "SummarizedExperiment"),
    function(x, modules, by = "rows", key = "row.names", as = "ids"){
    # Check margin
    if( length(by) == 0L || !by %in% c("rows", "cols") ){
        stop("'by' must be either 'rows' or 'cols'.", call. = FALSE)
    }
    # Select side information based on margin
    df <- if( by == "rows" ) rowData(x) else colData(x)
    # Build modules table
    out <- .get_modules(df, modules, key, as)
    # Convert to dataframe
    out <- as.data.frame(out)
    # Add tidy rownames
    rownames(out) <- rownames(df)
    return(out)
})


.get_modules <- function(df, modules, key, as){
    # Check if by is rownames
    is_rownames <- length(key) == 1L && key == "row.names"
    # Check args
    if( !is_rownames && !all(key %in% colnames(df)) ){
        stop("'key' must be 'row.names' or a variable of 'x'.", call. = FALSE)
    }
    if( !as %in% c("ids", "names") ){
        stop("'as' must be either 'ids' or 'names'.", call. = FALSE)
    }
    # Choose between ids and names
    target_col <- switch(as, ids = 2L, names = 3L)
    # Convert linkmap to wide format
    modules <- table(modules[c(1L, target_col)]) == 1
    # Match by rownames
    if( is_rownames ){
        idx <- match(rownames(df), rownames(modules))
    # Match by one variable
    }else if( length(key) == 1L ){
        idx <- match(df[[key]], rownames(modules))
    # Match by multiple variables (decreasing priority)
    }else{
        idx <- apply(df[key], 1L, function(row){
            m <- match(row, rownames(modules))
            first_match <- m[!is.na(m)][1L]
        })
    }
    # Select and order matched modules rows
    out <- modules[idx, , drop = FALSE]
    # Replace missing with false
    out[is.na(out)] <- FALSE
    return(out)
}
