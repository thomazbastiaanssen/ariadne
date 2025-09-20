#' Add information on modules of features or samples
#'
#' \code{getModules} and \code{addModules} generate a logical modules table from
#' a list of of named vectors specifying which features or samples belong to
#' each module.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param modules \code{Named character list}. Named list of vectors, where each
#'   vector is a module and its elements are the module members.
#'
#' @param by \code{Numeric scalar} or \code{Character scalar}. Determines if
#'   modules are added row-wise (\code{c(1L, "features")}) or column-wise.
#'   (\code{C(2L, "samples")}). (Default: \code{1L}).
#' 
#' @param group \code{Character vector}. One or more variable names from
#'   \code{rowData} or \code{colData}, depending on \code{by}. If set to
#'   \code{"taxonomy"}, modules are derived from the combined taxonomic ranks.
#'   (Default: \code{"taxonomy"})
#'
#' @param exact.tax.level \code{Logical scalar}. Should only the deepest
#'   taxonomic rank be used to determine whether a feature belongs to a module?
#'   if \code{FALSE}, all ranks are considered. (Default: \code{FALSE}).
#' 
#' @details
#' \code{getModules} and \code{addModules} allow to integrate information on
#' microbial modules into \code{rowData(x)} from databases like BugSigDB or
#' custom tables of microbes like \code{butyrate}. \code{sigs} must be a single
#' character vector or a list of character vectors, where each vector
#' corresponds to a module. Each microbial signature should follow the metaphlan
#' taxonomy format, e.g.,
#' \code{"k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Corynebacteriales"}.
#' 
#' @return
#' \code{getModules} returns a \code{matrix} where rows and columns represent
#' features/samples and modules, respectively. \code{addModules} returns an
#' updated \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object with the modules table in the \code{rowData} or \code{colData}.
#' 
#' @name getModules
#'
#' @examples
#' # Load butyrate module
#' data("butyrate", package = "ariadne")
#' 
#' # Load dataset
#' data("Tengeler2020", package = "mia")
#' tse <- Tengeler2020
#' 
#' # Convert butyrate table to module list
#' butyrate <- getFullTaxonomyLabels(butyrate)
#' 
#' # Generate modules table
#' mod.table <- getModules(tse, butyrate)
#' 
#' # Generate & store modules table
#' tse <- addModules(tse, modules)
#'
#' # Make example list of two modules
#' modules <- list(A = c("g__Bacteroides", "f__Enterobacteriaceae"),
#'                 B = c("c__Clostridia", "p__Cyanobacteria"))
#'
#' # Match modules with only the deepest taxonomic rank of tse
#' mod.table <- getModules(tse, modules, exact.tax.level = TRUE)
#' 
#' # Make example list of three column modules
#' col.modules <- list(A = c("Cohort_1", "Cohort_3"),
#'                     B = c("Cohort_2"),
#'                     c("Cohort_1", "Cohort_2", "Cohort_3"))
#' 
#' Find column modules based on cohort variable
#' mod.table <- getModules(tse, col.modules, by = 2L, group = "cohort")
#' 
NULL

#' @rdname getModules
#' @export
#' @importFrom SummarizedExperiment rowData colData
setMethod("addModules", signature = c(x = "SummarizedExperiment"),
    function(x, modules, by = 1L, group = "taxonomy", exact.tax.level = FALSE){
        # Make modules table
        modules <- getModules(
            x,
            modules,
            by = by,
            group = group,
            exact.tax.level = exact.tax.level
        )
        # Bind modules table with side information
        if( by %in% c(1L, "features") ){
            rowData(x) <- cbind(rowData(x), modules)
        }else if( by %in% c(2L, "samples") ){
            colData(x) <- cbind(colData(x), modules)
        }
        return(x)
    }
)

#' @export
#' @rdname getModules
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom mia taxonomyRanks
#' @importFrom stringr str_escape str_remove
setMethod("getModules", signature = c(x = "SummarizedExperiment"),
    function(x, modules, by = 1L, group = "taxonomy", exact.tax.level = FALSE){
        # Check modules
        if( !is.vector(modules) ){
            stop("'modules' must be a character vector or list of character ",
                "vectors, where each vector corresponds to a module.",
                call. = FALSE)
        }
        # Convert modules to list in case of only one module
        if( !is(modules, "list") ){
            modules <- list(module = modules)
        }
        # Check margin
        if( !by %in% c(1L, 2L, "features", "samples") ){
            stop("'by' must be one of 1, 2, features and samples.",
            call. = FALSE)
        }
        # Check exact.tax.level
        if( !is.logical(exact.tax.level) ){
            stop("'exact.tax.level' must be TRUE or FALSE.", call. = FALSE)
        }
        if( exact.tax.level && group != "taxonomy" ){
            warning("'exact.tax.label' is ignored when 'group' is not taxonomy.",
                call. = FALSE)
        }
        # Check group
        if( !is(x, "TreeSummarizedExperiment") && group == "taxonomy" ){
            stop("'group' can be 'taxonomy' only when 'x' is a ",
                "TreeSummarizedExperiment object.", call. = FALSE)
        }
        if( length(group) == 1L && group == "taxonomy" ){
            # Retrieve available taxrank colnames
            group <- taxonomyRanks(x)
            # Add rank prefixes to taxrank cols
            rowData(x)[ , group] <- .add_prefix_to_taxtable(rowData(x)[ , group])
            # Extract deepest taxonomic rank
            modules <- lapply(modules, str_remove, pattern = ".*\\|")
        }
        # Select side information
        x <- switch(by, rowData(x), colData(x))
        # Check group
        if( !all(group %in% names(x)) ){
            stop("'group' does not match any variables in the side information.",
                call. = FALSE)
        }
        # Collapse group variables to single strings
        group <- apply(x[ , group, drop = FALSE], 1L, function(row)
            paste(row, collapse = "|")
        )
        # Collapse signatures to single strings
        modules <- vapply(modules, function(mod){
            mod |>
                str_escape() |>
                paste0(collapse = "|")
            },
            character(1L)
        )
        # Reduce to deepest rank if exact.tax.level is on
        if( exact.tax.level ){
            group <- str_remove(group, ".*\\|")
        }
        # Make modules table
        mod.table <- .make_modules_table(group, modules)
        return(mod.table)
    }
)

# Define function to construct modules table based on bugsigdb signatures
#' @importFrom stringr str_detect
.make_modules_table <- function(group, modules){
    # Initialise all-FALSE modules table
    mod.table <- matrix(
        FALSE,
        nrow = length(group),
        ncol = length(modules)
    )
    # Find member features for each module
    for( i in seq_along(modules) ){
        mod.table[ , i] <- str_detect(group, modules[[i]])
    }
    # Add module names to table
    colnames(mod.table) <- names(modules)
    return(mod.table)
}

### HELPER FUNCTIONS ###

# Reduce taxcols of rowData to taxstring in metaphlan format
#' @export
#' @importFrom SummarizedExperiment rowData
getFullTaxonomyLabels <- function(x){
    # Add taxrank prefixes to taxcols of rowData
    tax <- .add_prefix_to_taxtable(x)
    # Collapse taxcols to taxstring in metaphlan format
    tax <- apply(tax, 1L, function(row) paste(row, collapse = "|"))
    # Remove empty taxranks
    tax <- gsub("(?:\\|[a-z]__)+$", "", tax)
    return(tax)
}

# Add taxrank prefixes to taxonomy table
.add_prefix_to_taxtable <- function(tax){
    # Retrieve taxrank prefixes
    tax.prefix <- paste0(mia:::getTaxonomyRankPrefixes()[tolower(names(tax))], "__")
    # Add prefix unless rank is NA
    tax <- mapply(function(rank, prefix)
        ifelse(is.na(rank), NA, paste0(prefix, rank)), tax, tax.prefix)
    return(tax)
}
