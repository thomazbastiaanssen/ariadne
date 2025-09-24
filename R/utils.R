#' Utility functions
#' 
#' These utility functions are used throughout the package and may be relevant
#' in other packages dealing with annotation mappings. \code{as.linkmap}
#' converts a list of named vectors to a linkmap data.frame.
#' 
#' @param tse A
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#'   object.
#' 
#' @param se A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @returns
#' \code{as.linkmap} returns a linkmap \code{data.frame} where the first and
#' second columns contains \code{keys} and \code{values} and each row represents
#' a unique combination of the two.
#'
#' @name utils
NULL

# Reduce taxcols of rowData to taxstring in metaphlan format
#' @export
#' @rdname utils
#' @importFrom SummarizedExperiment rowData
getFullTaxonomyLabels <- function(tse){
    # Add taxrank prefixes to taxcols of rowData
    tax <- .add_prefix_to_taxtable(tse)
    # Collapse taxcols to taxstring in metaphlan format
    tax <- apply(tax, 1L, paste, collapse = "|")
    # Remove empty taxranks
    tax <- gsub("(?:\\|[a-z]__)+$", "", tax)
    return(tax)
}

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
