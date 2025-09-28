#' @import S7
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL

#' @importFrom reticulate py_require import
.onLoad <- function(libname, pkgname){
    # Register S7 methods
    S7::methods_register()
    # Import Python dependencies
    reticulate::py_require("rdflib")
    rdflib <<- reticulate::import("rdflib")
    # Use unverified ssl for MacOS
    if( Sys.info()[["sysname"]] == "Darwin" ){
        reticulate::py_require("ssl")
        ssl <<- reticulate::import("ssl")
        ssl$`_create_default_https_context` <- ssl$`_create_unverified_context`
    }
}