#' @import S7
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL

rlang::on_load({
    S7::methods_register()
})

.onLoad <- function(...) {
    rlang::run_on_load()
}
