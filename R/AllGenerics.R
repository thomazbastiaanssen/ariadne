
igraph <- S7::new_S3_class("igraph")

#' @export
#' @rdname searchPath
searchPath <- S7::new_generic("searchPath", "graph")

#' @export
#' @rdname plotPath
plotPath <- S7::new_generic("plotPath", "graph")

#' @export
#' @rdname weavePath
weavePath <- S7::new_generic("weavePath", "graph")

#' @export
#' @rdname weavePath
weaveComplex <- S7::new_generic("weaveComplex", "graph")
