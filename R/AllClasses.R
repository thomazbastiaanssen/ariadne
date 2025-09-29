#' MultiFactor S7 container class
#' @name MultiFactorDB
#' @rdname MultiFactorDB-class
#' @description
#' `MultiFactorDB` is an S7 class to interface, organize and manage multiple
#' sets of factors in a remote database. `MultiFactorDB` behaves like a regular.
#' `MultiFactor`. Methods for `MultiFactorDB` aim to follow `factor` behaviour.
#'
#' @details
#' The most straightforward way to construct a `MultiFactor` object is as a
#' named list of named data.frames. The columns of the data.frames indicate the
#' category of factor in that column.
#'
#' A `MultiFactor` object presents itself similar to a `data.frame`, in the
#' sense that level types can be called as columns and individual data.frame
#' components can be called as rows.
#' `MultiFactor` inherits from `list`; Content can be accessed through regular
#' list methods (e.g., `[`, `[[`).
#' @slot levels `Named list of character vectors`. Accessed through `levels(x)`
#' @slot value Optional. `Named list of vectors` of same size as `MultiFactor`
#'     content, or `TRUE`, if missing.
#' @slot map `(sparse) Matrix` specifying which elements contain which levels.
#' @returns a `MultiFactor` object.
#' @export
#'
MultiFactorDB <- S7::new_class(
    "MultiFactorDB",
    package = "ariadne",
    parent = MultiFactor::MultiFactor,
    properties  = list(
        levels = S7::new_property(getter = function(self) lapply(self, levels)),
        value = S7::new_property(
            getter = function(self) lapply(self, \(x) x@value)
            ),
        map = S7::new_property(
            getter = function(self) MultiFactor:::.mapMultiFactor(self, mode = "pattern")
        )
    ),
    constructor = function(x) {

        S7::new_object(
            .parent = MultiFactor::MultiFactor(x)
        )
    },

)

#' LinkMapDB S7 container class
#' @name LinkMapDB
#' @rdname LinkMapDB-class
#' @description
#' `LinkMapDB`is an S7 class to interface, organize and manage sets of factors
#' in a remote database. Methods for `LinkMapDB` aim to follow `factor`
#' behaviour.
#'
#' @slot levels `Named list` of character vectors depicting levels.
#' @slot value `Character scalar` Repository information.
#' @slot is_bool `Boolean`, not used, for compatibility. Returns `TRUE`.
#' @param x a `data.frame` with two named columns that can be coerced to factors
#' @returns a `LinkMapDB` object.
#' @examples
#' LinkMapDB(ariadne:::ChocoPhlAn[[1L]], repo = "ChocoPhlAn")
#' @seealso [MultiFactorDB()]
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
    constructor = function(x, repo = c("ChocoPhlAn","Woltka")) {
        value <- match.arg(repo)
        S7::new_object(
            .parent = MultiFactor::LinkMap(x),
            value   = value,
            is_bool = TRUE
        )
    }
)

ChocoPhlAn = list(
    eggnog2uniref50    = data.frame(eggnog   = factor(), uniref50 = factor()),
    eggnog2uniref90    = data.frame(eggnog   = factor(), uniref90 = factor()),
    go2uniref50        = data.frame(go       = factor(), uniref50 = factor()),
    go2uniref90        = data.frame(go       = factor(), uniref90 = factor()),
    ko2uniref50        = data.frame(ko       = factor(), uniref50 = factor()),
    ko2uniref90        = data.frame(ko       = factor(), uniref90 = factor()),
    level4ec2uniref50  = data.frame(ec       = factor(), uniref50 = factor()),
    level4ec2uniref90  = data.frame(ec       = factor(), uniref90 = factor()),
    pfam2uniref50      = data.frame(pfam     = factor(), uniref50 = factor()),
    pfam2uniref90      = data.frame(pfam     = factor(), uniref90 = factor()),
    uniref502uniref90  = data.frame(uniref50 = factor(), uniref90 = factor())
)
