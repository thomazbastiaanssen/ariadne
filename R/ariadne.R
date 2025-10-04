#' Ariadne: On-demand annotation and module analysis
#' @name ariadne
#' @description
#' Display available prepared relational databases.
#' @importFrom rlang is_missing
#' @param repo `Missing` (default) or `Character or index scalar` selecting an
#'     available database
#' @returns if `repo` is missing, print names of available databases. Otherwise,
#'     return a `MultiFactor` object of the selected database.
#' @examples
#' ariadne()
#' ariadne("ChocoPhlAn")
#' @seealso [importMapping()]
#' @export
#'
ariadne <- function(repo) {
    repos <- c("ChocoPhlAn", "WoL")
    if(rlang::is_missing(repo)) {
        cat("Ariadne knows about the following relational databases:\n")
        cat({paste(repos, collapse = ", ")}, ".\n\n", sep = "")
        cat("Select them by name or index: `ariadne('WoL')` or `ariadne(1)`")
    }
    if(!rlang::is_missing(repo)) {
        x <- local({
            data("repo_layouts", package = "ariadne", envir = environment())
            return(mget(repos))
        })[[repo]]
        MultiFactor::MultiFactor(lapply(x, LinkMapDB))
    }
}
