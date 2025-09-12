#' @importFrom reticulate py_require import
.onLoad <- function(libname, pkgname){
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