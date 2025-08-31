#' @importFrom reticulate py_require
.onLoad <- function(libname, pkgname){
    py_require("rdflib")
    py_require("ssl")
}