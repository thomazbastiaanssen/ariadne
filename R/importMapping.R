#' Import mappings from a file or a database
#' 
#' \code{importMapping} retrieves mapping information from a file or a database.
#' 
#' @param map.file \code{Character vector}. One or more paths to custom mapping
#'   files or one or more of the available databases (\code{"ChocoPhlAn"}).
#'
#' @param from \code{Character vector}. One or more strings specifying the
#'   keys from which mapping is performed. (Default: \code{NULL})
#'
#' @param to \code{Character vector}. One or more strings specifying the
#'   values to which mapping is performed. (Default: \code{NULL})
#'
#' @param merge \code{Logical scalar}. Should multiple mapping files be merged.
#'   (Default: \code{TRUE}).
#' 
#' @param message \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#' 
#' @details
#' Structure of \code{map.file}
#' 
#' Currently, the following databases are available:
#' \itemize{
#'   \item{\code{"ChocoPhlAn"}: x-to-uniref mapping files from bioBakery
#'     \itemize{
#'       \item{repo: \href{https://zenodo.org/records/17100034}{https://zenodo.org/records/17100034}}
#'       \item{from: \code{c("eggnog", "go", "ko", "level4ec")}}
#'       \item{to: \code{c("uniref50", "uniref90")}}}
#'   }
#'   \item{\code{"WoLtka"}: description
#'     \itemize{
#'       \item{repo: \href{https://ftp.microbio.me/pub/wol-20April2021/}{https://ftp.microbio.me/pub/wol-20April2021/}}
#'       \item{from: }
#'       \item{to: }
#'     }}
#' }
#' 
#' @return
#' \code{importMapping} returns a named list of vectors, where each vector is a
#' mapping key and its elements are the mapped values.
#'
#' @examples
#' # Import eggnog-to-uniref90 mapping from ChocoPhlAn
#' map1 <- importMapping("ChocoPhlAn", from = "eggnog", to = "uniref90")
#' 
#' # Import ko-to-uniref90 mapping from ChocoPhlAn
#' map2 <- importMapping("ChocoPhlAn", from = "ko", to = "uniref90")
#' 
#' # Import and merge ko- and eggnog-to-uniref90 mappings from ChocoPhlAn
#' map3 <- importMapping(
#'     "ChocoPhlAn",
#'     from = c("eggnog", "ko"),
#'     to = "uniref90"
#' )
#' 
#' # Import local mapping file
#' # map4 <- importMapping("path/to/file")
#' 
#' @name importMapping
NULL

MappingDatabases <- list(
    ChocoPhlAn = list(
        repo = "https://zenodo.org/records/17100034/files/",
        from = c("eggnog", "go", "ko", "level4ec"),
        to = c("uniref50", "uniref90"),
        prefix = "map_",
        suffix = ".txt.gz",
        sep = "_"
    )
    # WoLtka <- list()
)

#' @rdname importMapping
#' @export
setMethod("importMapping", signature = c(map.file = "character"),
    function(map.file, from = NULL, to = NULL, merge = TRUE, message = TRUE){
        # Find number of mapping files
        input.lengths <- lengths(list(map.file, from, to))
        max.length <- max(input.lengths)
        # Check arguments
        if( !all(input.lengths %in% c(1, max.length)) ){
            stop("'map.files', 'from' and 'to' must have compatible lengths.",
                call. = FALSE)
        }
        if( !is.logical(merge) ){
            stop("'merge' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !is.logical(message) ){
            stop("'message' must be TRUE or FALSE.", call. = FALSE)
        }
        # Import each mapping file
        map <- mapply(
            .import_mapping,
            map.file = map.file, from = from, to = to,
            MoreArgs = list(message = message),
            SIMPLIFY = FALSE, USE.NAMES = FALSE
        )
        # Merge mapping files
        if( merge ){
            map <- unlist(map, recursive = FALSE)
        }
        return(map)
    }
)

# Import single mapping file
.import_mapping <- function(map.file, from, to, message){
    # Whether to use package or custom mapping
    if( map.file %in% names(MappingDatabases) ){
        # Construct path to database
        map.file <- .getPath(map.file, from, to)
        # Cache database
        map.file <- .getCache(map.file)
    }
    if( message ){
        message("Retrieving mappings from ", map.file, ".")
    }
    # Read file content
    lines <- readLines(map.file)
    # Split elements in each line by tab
    items <- strsplit(x = lines, split = "\t")
    # Extract keys
    keys <- vapply(items, FUN = function(x) x[1], FUN.VALUE = "")
    # Extract values
    values <- lapply(items, FUN = function(x) x[-1])
    # Add names to list
    names(values) <- keys
    return(values)
}

# Retrieve mapping file from database
.getPath <- function(map.file, from, to){
    # Retrieve database
    map.db <- MappingDatabases[[map.file]]
    if( !is.null(from) && !from %in% map.db$from ){
        stop("'from' should be defined and be one of ",
            paste(map.db$from, collapse = ", "), ".", call. = FALSE)
    }
    if( !is.null(to) && !to %in% map.db$to ){
        stop("'to' should be defined and be one of ",
            paste(map.db$to, collapse = ", "), ".", call. = FALSE)
    }
    map.file <- paste0(
        map.db$repo, map.db$prefix, from, map.db$sep, to, map.db$suffix
    )
    return(map.file)
}
