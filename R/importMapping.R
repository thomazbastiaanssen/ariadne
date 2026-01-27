#' Import mappings from a file or a database
#'
#' @name importMapping
#' @rdname importMapping
#'
#' @description
#' \code{importMapping} retrieves mapping information from a file or a database.
#'
#' @param x \code{Character vector}. One or more paths to custom mapping
#'   files or one or more names from the available databases
#'   (\code{c("ChocoPhlAn", "Woltka")}).
#'
#' @param subset \code{formula scalar} or \code{Character vector}. Specifies the
#'   'to' and 'from' of mapping in the syntax \code{from ~ to}, or equivalently
#'   \code{c("from", "to")}. (Default: \code{NULL})
#' @param verbose \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#'
#' @details
#'
#' Currently, the following databases are available:
#' \itemize{
#'   \item{\code{"ChocoPhlAn"}: x-to-uniref mapping files from bioBakery
#'     \itemize{
#'       \item{repo: \href{https://zenodo.org/records/17100034}{https://zenodo.org/records/17100034}}
#'       \item{from: \code{c("eggnog", "go", "ko", "level4ec")}}
#'       \item{to: \code{c("uniref50", "uniref90")}}}
#'   }
#'   \item{\code{"Woltka"}: x-to-uniref mapping files from the Web of Life
#'     \itemize{
#'       \item{repo: \href{https://ftp.microbio.me/pub/wol-20April2021/}{https://ftp.microbio.me/pub/wol-20April2021/}}
#'       \item{from: c("eggnog", "go", "ko", "orthodb", "refseq")}
#'       \item{to: "uniref90"}
#'     }}
#' }
#'
#' @return
#' \code{importMapping} returns a named list of vectors, where each vector is a
#' mapping key and its elements are the mapped values.
#'
#' @examples
#' # Import eggnog-to-uniref90 mapping from ChocoPhlAn
#' map1 <- importMapping("ChocoPhlAn", eggnog ~ uniref90)
#'
#' # Import ko-to-uniref90 mapping from ChocoPhlAn
#' map2 <- importMapping("ChocoPhlAn", ko ~ uniref90)
#'
#' # Import several files from an ariadne-indexed database:
#' db <- ariadne("ChocoPhlAn", ko ~ uniref90)
#'
#' # Which db-side files would be downloaded? Download them if not dry.run.
#' db.local <- importMapping(db, dry.run = TRUE)
#'
NULL

#' @importFrom MultiFactor LinkMap MultiFactor as.LinkMap
S7::method(importMapping, MultiFactor) <-
    function( x, subset = NULL, dry.run = TRUE ) {

    x <- subset(x, subset)
    lmdbs <-
        names(x)[vapply( x, inherits, "ariadne::LinkMapDB", FUN.VALUE = FALSE )]

    if(length(lmdbs) == 0L) {
        cat("This MultiFactor is up-to-date. No need to import anything. ")
        return(invisible(x))
    }
    if(dry.run) {
        cat("Disabling `dry.run` would download the following linkage files:\n")
        cat(paste(lmdbs, collapse = ", "))
    } else {
        cat("Downloading the following linkage files:\n")
        cat(paste(lmdbs, collapse = ", "))
        MultiFactor::MultiFactor(lapply(x, as.LinkMap))
    }
}

#' @importFrom MultiFactor LinkMap MultiFactor as.LinkMap
S7::method(importMapping, S7::class_character) <-
    function( x, subset = NULL, dry.run = TRUE ){

    importMapping(
        x = ariadne(x),
        subset = subset,
        dry.run = dry.run
    )
}

# Import single mapping file
.import_mapping <- function(x, verbose = TRUE){
    x <- .getCache(x)
    if( verbose ){
        message("Retrieving mappings from ", x, ".")
    }

    # Read file content
    line.content <- readLines(x)

    # Split elements in each line by tab
    line.content <- strsplit(line.content, "\t", fixed = TRUE)
    # Extract keys
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = character(1L))
    # Extract values
    values <- lapply(line.content, `[`, -1L)

    data.frame(
        id.x = rep(keys, lengths(values, use.names = FALSE)),
        id.y = unlist(values, recursive = TRUE, use.names = FALSE)
    )
}


#' @importFrom MultiFactor as.LinkMap
.process_woltka <- function(woltka.map){
    names(woltka.map) <- paste0("UniRef90_", names(woltka.map))
    linkmap <- as.LinkMap(woltka.map)
    woltka.map <- split(linkmap$x, linkmap$y)
    return(woltka.map)
}


# Read a df, make a MF-shaped list
.reftableToDFList <- function(x) `names<-`(lapply(
    seq_len(NROW(x)),
    FUN = function(y) `names<-`(
        data.frame(factor(),
                   factor()),
        c(x[y, 1:2])
    )),     paste(x[[1L]], x[[2L]], sep = "2")

)

# Read a MF-shaped list, make a df.
.ListToReftable <- function(x) `colnames<-`(
    as.data.frame(t(vapply(x, names, c("", "")))), c("from", "to")
)
