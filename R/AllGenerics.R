
#' @export
#' @rdname importMapping
importMapping <- S7::new_generic("importMapping", "map.file")

#' @export
#' @rdname importModules
importModules <- S7::new_generic("importModules", "module.file")

#' @export
#' @rdname mapModules
mapModules <- S7::new_generic("mapModules", "modules")

#' @export
#' @rdname getModules
getModules <- S7::new_generic("getModules", "x")

#' @export
#' @rdname getModules
addModules <- S7::new_generic("addModules", "x")
