#' @title Make into a LinkMap
#' @name as.LinkMap.LinkMapDB
#' @description
#' Return a LinkMap with the properties advertised in the `LinkMapDB`.
#' @return a `LinkMap` object
#' @importFrom MultiFactor as.LinkMap
#' @examples
#' # generate a MultiFactor
#' MF.db <- MultiFactor::MultiFactor(ariadne("ChocoPhlAn"))
#'
#' # Extract one underlying LinkMapDB object
#' LM.db <- MF.db[["ko2uniref90"]]
#'
#' LM <- MultiFactor::as.LinkMap(LM.db)
NULL

#' @export
#'
S7::method(as.LinkMap, LinkMapDB) <- function(x) {

    outnames <- names(x)

    x <- .urlLinkMapDB(x)
    x <- .import_mapping(x)

    names(x) <- outnames

    MultiFactor::LinkMap(x)
}


#' @param x LinkMapDB
#' @noRd
.urlLinkMapDB <- function(x) {
    switch (
        x@value,
        ChocoPhlAn = paste0(
            "https://zenodo.org/records/17100034/files/map_",
            names(x)[1L], "_", names(x)[2L], ".txt.gz"
        ),
        WoL = paste0(
            "https://ftp.microbio.me/pub/wol-20April2021/function",
            names(x)[1L], "/", names(x)[2L], ".map.xz"
        ),
        GO = paste0(
            "https://current.geneontology.org/ontology/external2go/",
            names(x)[1L], "2", names(x)[2L]
        )
    )
}
