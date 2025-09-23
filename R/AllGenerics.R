#' @rdname importMapping
#' @export
setGeneric("importMapping", signature = "map.file",
    function(map.file, from = NULL, to = NULL, merge = TRUE, verbose = TRUE)
    standardGeneric("importMapping"))

#' @rdname importModules
#' @export
setGeneric("importModules", signature = "module.file",
    function(module.file, mode = "single", merge = TRUE, verbose = TRUE)
    standardGeneric("importModules"))

#' @rdname mapModules
#' @export
setGeneric("mapModules", signature = "modules",
    function(modules, map, mode = "single", remove.empty = TRUE, verbose = TRUE)
    standardGeneric("mapModules"))

#' @rdname queryTaxonomy
#' @export
setGeneric("queryTaxonomy", signature = "map",
    function(map, remove.empty = TRUE, verbose = TRUE)
    standardGeneric("queryTaxonomy"))

#' @rdname getModules
#' @export
setGeneric("getModules", signature = "x",
    function(x, modules, by = 1L, group = "taxonomy", exact.tax.level = FALSE)
    standardGeneric("getModules"))

#' @rdname getModules
#' @export
setGeneric("addModules", signature = "x",
    function(x, modules, by = 1L, group = "taxonomy", exact.tax.level = FALSE)
    standardGeneric("addModules"))

#' @rdname utils
#' @export
setGeneric("as.linkmap", signature = "values",
    function(values, keys = NULL, col.names = NULL)
    standardGeneric("as.linkmap"))
