#' List resource versions registered in the ariadne database
#' 
#' @name listResourceVersions
#' 
#' @description
#' listResourceVersions shows the available versions for the different resources
#' registered in the ariadne database. Any of those versions can be passed to
#' \code{\link{ariadne}}. Resources labelled with \code{"latest"} cannot be
#' versioned as they are dynamically accessed via their API or SPARQL endpoint.
#' 
#' @param default \code{Logical scalar}. Should only default versions be listed.
#'   (Default: \code{FALSE})
#' 
#' @returns A data.frame with information on registered resource versions.
#' 
#' @seealso \code{\link{ariadne}}
#' 
#' @examples
#' # View all available resource versions
#' listResourceVersions()
#' 
#' # View default resource versions
#' listResourceVersions(default = TRUE)
NULL

#' @export
#' @rdname listResourceVersions
listResourceVersions <- function(default = FALSE){
    # Retrieve metadata on resource versions
    meta <- versionMetadata
    # If default is turned on
    if( default ){
        # Select only default versions
        meta <- meta[meta$default, ]
    }
    # Build resource base urls
    urls <- attr(meta, "urls")
    meta$url <- urls[match(meta$source, names(urls))]
    meta$url <- paste0(meta$url, meta$key, "/")
    # Rename source column to resource
    meta$resource <- meta$source
    # Select relevant columns to print
    meta <- meta[ , c("resource", "version", "url")]
    return(meta)
}


#' Process HUMAnN gene families
#' 
#' @name processGeneFamilies
#' 
#' @description
#' processGeneFamilies prepares a SummarizedExperiment object containing the
#' HUMAnN gene families, such as those provided by curatedMetagenomicData, so
#' its feature-wise information on genes and taxa are added to the rowData. This
#' makes ariadne interoperable with HUMAnN gene families data.
#' 
#' @param x A
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#' 
#' @returns
#' An object of the same type as \code{x} with four additional columns in the
#' rowData: uniref90, taxname, genus and species.
#' 
#' @examples
#' library(curatedMetagenomicData)
#' 
#' # Import gene families
#' genes <- curatedMetagenomicData(
#'     "AsnicarF_2017.gene_families",
#'     dryrun = FALSE
#' )
#' 
#' # Extract experiment from list
#' genes <- genes[[1]]
#' 
#' # Process gene families
#' genes <- processGeneFamilies(genes)
#' 
#' # Print head of rowData
#' head(rowData(genes, use.names = FALSE))
NULL

#' @export
#' @rdname processGeneFamilies
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stringr fixed str_detect str_split
processGeneFamilies <- function(x){
    # Select rows with non-null taxa
    x <- x[str_detect(rownames(x), fixed("|")), ]
    x <- x[str_detect(rownames(x), "unclassified", negate = TRUE), ]
    # Split gene and taxonomy
    gene.linkmap <- as.data.frame(
        str_split(rownames(x), fixed("|"), n = 2, simplify = TRUE)
    )
    names(gene.linkmap) <- c("uniref90", "taxname")
    # Split genus and species
    tax.linkmap <- as.data.frame(
        str_split(gene.linkmap$taxname, fixed("."), n = 2, simplify = TRUE),
    )
    names(tax.linkmap) <- c("genus", "species")
    # Bind gene and tax linkmaps
    rowData(x) <- cbind(rowData(x), gene.linkmap, tax.linkmap)
    return(x)
}
