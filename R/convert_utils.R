#' Convert a MultiFactor to three-column data.frame format.
#' @param x MultiFactor
#' @param source_name `Character scalar` Content of source column in output.
#' @returns a `data.frame` with three columns; to, from and source.
#' @examples
#' data("ChocoPhlAn_layout",  package = "ariadne")
#' data("WoL_layout",  package = "ariadne")
#'
#' # Three-column dfs of layouts
#' cho_df <- MF2DF(ChocoPhlAn)
#' wol_df <- MF2DF(WoL)
#'
#' # Stack them into tall three-column df.
#' tot_df <- rbind(cho_df, wol_df)
#'
#' DF2MF(tot_df)
#'
MF2DF <- function(x, source_name) {
    if(missing(source_name)) source_name <- deparse1(substitute(x))
    x <- MultiFactor::MultiFactor(x)
    x <- lapply(x, names) |>
        list2DF() |>
        t() |>
        as.data.frame.matrix(row.names = NULL, make.names = FALSE)
    rownames(x) <- NULL
    colnames(x) <- c("from", "to")
    x[["source_name"]] <- source_name
    x
}

#' Build a MultiFactor from three-columns data.frame notation.
#' @param x `data.frame` with three columns; to, from and source.
#' @returns a MultiFactor.
#' @rdname MF2DF
#'
DF2MF <- function(x) {
    mf <- tapply(x, x[[3L]], .snipToLinkMap, simplify = FALSE)

    MultiFactor::MultiFactor(unlist(mf, recursive = FALSE))
}

.snipToLinkMap <- function(x) {
    out <- apply(
        x, MARGIN = 1L,
        function(x) LinkMap( `colnames<-`(
            data.frame("<placeholder>","<placeholder>"), x[seq_len(2L)]
        ),
        metadata = list(source = x[3L])
        ), simplify = FALSE
    )
    out
}


