#' @title Make into a LinkMap
#' @name as.LinkMap.LinkMapDB
#' @description
#' Return a LinkMap with the properties advertised in the `LinkMapDB`.
#' @return a `LinkMap` object
#' @importFrom MultiFactor as.LinkMap
#' @examples
#' # generate a MultiFactorDB
#' MF.db <- MultiFactorDB(ariadne:::ChocoPhlAn)
#'
#' # Extract one underlying LinkMapDB object
#' LM.db <- MF.db[["ko2uniref90"]]
#'
#' LM <- MultiFactor::as.LinkMap(LM.db)
NULL

#' @export
#'
S7::`method<-`(MultiFactor::as.LinkMap, LinkMapDB, function(x) {

    outnames <- names(x)
    x <- .getCache(.urlLinkMapDB(x))

    # Read file content
    line.content <- readLines(x)

    # Split elements in each line by tab
    line.content <- strsplit(line.content, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, FUN = function(x) x[1], FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line.content, FUN = function(x) x[-1])

    x <- `names<-`(
        data.frame(
            id.x = rep(keys, vapply(values, length, 1L)),
            id.y = unlist(values, recursive = TRUE, use.names = FALSE)
        ), outnames
    )
    MultiFactor::LinkMap(x)
}
)


#' @param x LinkMapDB
#' @noRd
.urlLinkMapDB <- function(x) {
    switch (
        x@value,
        ChocoPhlAn = paste0(
            "https://zenodo.org/records/17100034/files/map_",
            names(x)[1L], "_", names(x)[2L], ".txt.gz"
        ),
        Woltka     = paste0(
            "https://ftp.microbio.me/pub/wol-20April2021/function",
            names(x)[1L], "/", names(x)[2L], ".map.xz"
        )
    )
}
