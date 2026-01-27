#' Ariadne: On-demand annotation and module analysis
#' @name ariadne
#' @description
#' Display available prepared relational databases.
#' @importFrom rlang is_missing
#' @importFrom rlang dots_list
#' @param x `Missing` (default) or `Character vector` selecting any number of
#'     available relational database or module set.
#' @param .by Either a `formula` or a `Character vector` of length 2 with the
#'     names of the desired combination of feature types. See `MultiFactor`
#'     subset method, which is called internally.
#' @returns if `x` is missing, print names of available databases. Otherwise,
#'     return a `MultiFactor` object of the selected database.
#' @examples
#' ariadne()
#' ariadne("Choco")
#' ariadne("ChocoPhlAn", ko ~ uniref90)
#' ariadne(c("ChocoPhlAn", "GBM"), GBM ~ uniref90)
#' @seealso [importMapping()]
#' @seealso [importModules()]
#' @export
#'
ariadne <- function(x, .by) {
    repos   <- c("ChocoPhlAn", "WoL")
    modules <- c("GBM", "GMM")

    if(rlang::is_missing(x)) {
        stopifnot(
            "'.by' can't be used without specifying 'x'. " =
                      rlang::is_missing(.by)
            )
        return(.ariadne_report(repos, modules))
    }

    res <- .ariadne_single(repos, modules, x)

    if(!rlang::is_missing(.by)) {
        res <- subset(res, .by)
    }
    return(res)
}


.ariadne_report <-function(repos, modules) {
    cat("Ariadne knows about the following relational databases:\n")
    cat("  ", {paste(repos, collapse = ", ")}, "\n", sep = "")
    cat("and the following module sets:\n")
    cat("  ", {paste(modules, collapse = ", ")}, ".\n\n", sep = "")
    cat("Select them by name: `ariadne('WoL')`")
    return(invisible())
}


.ariadne_single <- function(repos, modules, x) {

    x <- c(repos, modules)[pmatch(
        toupper(x),
        toupper(c(repos, modules)),
        duplicates.ok = FALSE
    )]
    
    if(all(is.na(x))) {
        stop("'x' should select a known relational database or module set:\n",
            "  ", {paste(c(repos, modules), collapse = ", ")}, ".",
            call. = FALSE)
    }
    
    do.call(`c`, lapply(x, function(x) {
        if ( x %in% repos ) return(.get_repo(x))
        if ( x %in% modules ) return(importModules(x))
    }))
}
