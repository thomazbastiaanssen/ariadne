
setOldClass("igraph")

#' @export
#' @rdname searchPath
setGeneric("searchPath", signature = c("graph"), function(graph, ...)
    standardGeneric("searchPath"))

#' @export
#' @rdname plotPath
setGeneric("plotPath", signature = c("graph"), function(graph, ...)
    standardGeneric("plotPath"))

#' @export
#' @rdname drawPath
setGeneric("drawPath", signature = c("graph"), function(graph, ...)
    standardGeneric("drawPath"))

#' @export
#' @rdname linkNames
setGeneric("linkNames", signature = c("graph"), function(graph, ...)
    standardGeneric("linkNames"))

#' @export
#' @rdname addResource
setGeneric("addResource", signature = c("graph"), function(graph, ...)
    standardGeneric("addResource"))

#' @export
#' @rdname weavePath
setGeneric("weavePath", signature = c("graph"), function(graph, ...)
    standardGeneric("weavePath"))

#' @export
#' @rdname weavePath
setGeneric("weaveComplex", signature = c("graph"), function(graph, ...)
    standardGeneric("weaveComplex"))

#' @export
#' @rdname recallPath
setGeneric("recallPath", signature = c("path_df"), function(path_df, ...)
    standardGeneric("recallPath"))

#' @export
#' @rdname addModules
setGeneric("getModules", signature = c("x"), function(x, ...)
    standardGeneric("getModules"))

#' @export
#' @rdname addModules
setGeneric("addModules", signature = c("x"), function(x, ...)
    standardGeneric("addModules"))
