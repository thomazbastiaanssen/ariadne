#' Ariadne: On-demand annotation and module analysis
#' @name ariadne
#' @description
#' Display available prepared relational databases.
#' @importFrom rlang is_missing
#' @param x `Missing` (default) or `Character or index scalar` selecting an
#'     available relational database or module set.
#' @returns if `x` is missing, print names of available databases. Otherwise,
#'     return a `MultiFactor` object of the selected database.
#' @examples
#' ariadne()
#' ariadne("ChocoPhlAn")
#' @seealso [importMapping()]
#' @seealso [importModules()]
#' @export
#'
ariadne <- function(x) {
    repos   <- c("ChocoPhlAn", "WoL")
    modules <- c("GBM",        "GMM")

    if(rlang::is_missing(x)) {
        cat("Ariadne knows about the following relational databases:\n")
        cat("  ", {paste(repos, collapse = ", ")}, "\n", sep = "")
        cat("and the following module sets:\n")
        cat("  ", {paste(modules, collapse = ", ")}, ".\n\n", sep = "")
        cat("Select them by name or index: `ariadne('WoL')` or `ariadne(1)`")
        return(invisible())
    }

    if(is.numeric(x)) x <- c(repos, modules)[x]
    x <- c(repos, modules)[pmatch(
        toupper(x),
        toupper(c(repos, modules)),
        duplicates.ok = FALSE
    )]
    if(!x %in% c(repos, modules)) {
        stop(
            "'x' should select a known relational database or module set:\n",
            "  ", {paste(c(repos, modules), collapse = ", ")}, ".")
    }
    if(x %in% repos) {
        out <- local({
            data("repo_layouts", package = "ariadne", envir = environment())
            return(get(x))
        })

        out <- MultiFactor::MultiFactor(lapply(out, LinkMapDB))
        return(out)
    }
    if (x %in% modules) {
        out <- importModules(x)
        return(out)
    }

}
