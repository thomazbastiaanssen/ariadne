#' @rdname importMapping
#' @export
setGeneric("importMapping", signature = "map.file",
    function(map.file, from = NULL, to = NULL, merge = TRUE, message = TRUE)
    standardGeneric("importMapping"))

#' @rdname importModules
#' @export
setGeneric("importModules", signature = "module.file",
    function(module.file, merge = TRUE, message = TRUE)
    standardGeneric("importModules"))

#' @rdname mapModules
#' @export
setGeneric("mapModules", signature = "modules",
    function(modules, map, remove.empty = TRUE, message = TRUE)
    standardGeneric("mapModules"))

#' @rdname utils
#' @export
setGeneric("as.linkmap", signature = "values",
    function(values, keys = NULL, col.names = NULL)
    standardGeneric("as.linkmap"))
