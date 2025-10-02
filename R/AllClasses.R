#' LinkMapDB S7 container class
#' @name LinkMapDB
#' @rdname LinkMapDB-class
#' @description
#' `LinkMapDB`is an S7 class to interface, organize and manage sets of factors
#' in a remote database. Methods for `LinkMapDB` aim to follow `factor`
#' behaviour. The user is not expected to interface with `LinkMapDB` objects
#' directly.
#'
#' @slot levels `Named list` of character vectors depicting levels.
#' @slot value `Character scalar` Repository information.
#' @slot is_bool `Boolean`, not used, for compatibility. Returns `TRUE`.
#' @param x a `data.frame` with two named columns that can be coerced to factors
#' @returns a `LinkMapDB` object.
#' @examples
#' LinkMapDB
#' @export
#'
LinkMapDB <- S7::new_class(
    "LinkMapDB",
    package = "ariadne",
    parent = MultiFactor::LinkMap,
    properties  = list(
        levels   = S7::new_property(getter = function(self) lapply(self, levels)),
        value    = S7::class_character,
        is_bool  = S7::class_logical
    ),
    constructor = function(x, repo = c("ChocoPhlAn","WoL", "other")) {
        value <- match.arg(repo)
        S7::new_object(
            .parent = MultiFactor::LinkMap(x),
            value   = value,
            is_bool = TRUE
        )
    }
)

