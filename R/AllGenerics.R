
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
#' @rdname weavePath
setGeneric("weavePath", signature = c("graph"), function(graph, ...)
    standardGeneric("weavePath"))

#' @export
#' @rdname weavePath
setGeneric("weaveComplex", signature = c("graph"), function(graph, ...)
    standardGeneric("weaveComplex"))
