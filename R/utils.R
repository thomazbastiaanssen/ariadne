#' List resource versions registered in the ariadne database
#' 
#' @name listResourceVersions
#' @rdname listResourceVersions
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
#' @examples
#' # View all available resource versions
#' listResourceVersions()
#' 
#' # View default resource versions
#' listResourceVersions(default = TRUE)
#' 
#' @seealso \code{\link{ariadne}}
#' 
#' @export
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
#' @rdname processGeneFamilies
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
#' 
#' @export
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


#' Append linkmaps to SummarizedExperiment side information
#' 
#' @name appendModules
#' @rdname appendModules
#' 
#' @description
#' appendModules allows to add the output of weavePath or complexPath to the
#' rowData or colData of a SummarizedExperiment (SE) object. This makes ariadne
#' interoperable with SE-based data analysis.
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
#' @param as \code{Character scalar} A string specifying whether target ids or
#'   names should be appended to \code{x}. (Default: \code{"ids"})
#' 
#' @returns
#' An object of the same type as \code{x} with additional columns, each
#' containing information on membership of a feature to a certain module.
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
#' # Add butyrate-producer module to rowData
#' rowData(tse) <- appendModules(rowData(tse), butyrate, by = "Genus")
#' 
#' # Add modules based on multiple variables given in order of priority
#' rowData(tse) <- appendModules(
#'     rowData(tse), butyrate, by = c("Genus", "Family")
#' )
#' 
#' # Add module names instead of ids
#' rowData(tse) <- appendModules(
#'     rowData(tse), butyrate, by = "Genus", as = "names"
#' )
#' 
#' # Generate relative abundance table
#' tse <- transformAssay(tse, method = "relabundance")
#' 
#' # Agglomerate features by membership to butyrate-producer module
#' modules <- agglomerateByModule(tse, by = "rows", group = "butyrate")
#' 
#' # Plot relative abundance of butyrate producers
#' plotAbundance(modules, assay.type = "relabundance")
#' 
#' @export
appendModules <- function(x, modules, by = "row.names", as = "ids"){
    # Check if by is rownames
    is_rownames <- length(by) == 1L && by == "row.names"
    # Check args
    if( !is_rownames && !all(by %in% colnames(x)) ){
        stop("'by' must be 'row.names' or a variable of 'x'.", call. = FALSE)
    }
    if( !as %in% c("ids", "names") ){
        stop("'as' must be either 'ids' or 'names'.", call. = FALSE)
    }
    # Choose between ids and names
    target_col <- switch(as, ids = 2L, names = 3L)
    # Convert linkmap to wide format
    modules <- table(modules[c(1L, target_col)]) == 1
    
    if( is_rownames ){
        idx <- match(rownames(x), rownames(modules))
    }else if( length(by) == 1L ){
        idx <- match(x[[by]], rownames(modules))
    }else{
        idx <- apply(x[by], 1L, function(row){
            m <- match(row, rownames(modules))
            first_match <- m[!is.na(m)][1L]
        })
    }
    
    modules <- modules[idx, , drop = FALSE]
    modules[is.na(modules)] <- FALSE
    
    to_remove <- which(colnames(x) %in% colnames(modules))
    
    if( length(to_remove) != 0L ){
        warning("Some columns of 'x' were replaced.", call. = FALSE)
        x[to_remove] <- NULL
    }
    
    out <- cbind(x, modules)
    return(out)
}
