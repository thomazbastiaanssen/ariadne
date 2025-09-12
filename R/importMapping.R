#' Import mappings from a file or a database
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
    function(map.file, from = NULL, to = NULL, message = TRUE){
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
)

#' @rdname importMapping
#' @export
setMethod("importMappings", signature = c(map.files = "character"),
    function(map.files, from = NULL, to = NULL, merge = TRUE, message = TRUE){
        if( !is.null(from) && length(map.files) != length(from) ){
            stop("'from' should have the same length as 'map.files")
        }
        if( !is.null(to) && length(map.files) != length(to) ){
            stop("'to' should have the same length as 'map.files")
        }
        if( !is.logical(merge) ){
            stop("'merge' must be TRUE or FALSE.", call. = FALSE)
        }
        # Import each mapping file
        maps <- mapply(
            importMapping,
            map.files = map.files, from = from, to = to,
            MoreArgs = list(message = message)
        )
        # Merge mapping files
        if( merge ){
            maps <- unlist(maps, recursive = FALSE)
        }
        return(maps)
    }
)
    
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
