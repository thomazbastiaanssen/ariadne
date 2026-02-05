#' Map modules to taxa
#' @name mapModules
#' @rdname mapModules
#'
#' @description
#' \code{mapModules} returns a list of modules containing the taxa that are
#' members of each module. Taxa are derived from uniref ids by querying UniProt
#' SPARQL. Membership is based on whether a taxon meets the criteria specified
#' by a module.
#'
#' @param map \code{MultiFactor}. Typically produced by `importModules()`.
#'
#' @param x \code{matrix}. Feature table that can be coerced to `matrix`, such
#' as a `data.frame`, with features as rows and samples as columns.
#'
#' @param x.type \code{Character scalar}. Type of ID in \code{row.names(x)}.
#'     (Default: \code{"uniref90"}).
#'
#' @param x.stratfied \code{Logical or character scalar}. If \code{x} is NOT a
#'     stratified table; `FALSE` (Default). Otherwise, if \code{x} is stratified
#'     by rows, provide the field separator character.
#'     For example: \code{feature|subtype} requires \code{x.stratified = "|"}.
#' @param method \code{Character vector}. Define what output parameter should
#'     be computed. One of \code{c("sum", "count", "coverage", "presence")}.
#'     (Default: \code{"sum"}).
#'
#' @param coverage.threshold \code{Numeric scalar}. Minimum proportion of
#'     components (i.e., reaction steps) required to be considered present.
#'     bounded between 0-1. (Default: 0.8)
#'
#' @param verbose \code{Logical scalar}. Should information on execution be
#'     printed in the console. (Default: \code{TRUE}).
#'
#' @returns \code{mapModules} returns a numeric matrix with modules as rows and
#'     samples as columns, with content depending on desired `mode` argument.
#'
#' @examples
#' # Import GBM
#' map <- importModules("GBM")
#'
#' # sparse feature table with 20 random samples
#' x_present <- replicate(20,
#'     sample(
#'     c(TRUE, FALSE), size = nlevels(map$ko_complex2ko[[2L]]), replace = TRUE
#'     ),
#'     simplify = TRUE
#' )
#' x <- x_present * matrix(
#'     rnorm(n = prod(dim(x_present))),
#'     nrow = NROW(x_present), ncol = 20
#' )^2
#'
#' # Map modules to feature table
#' mapModules(map, x, x.type = "ko", method = "coverage")
#' mapModules(map, x, x.type = "ko", method = "sum")
#' mapModules(map, x, x.type = "ko", method = "count")
#' mapModules(map, x, x.type = "ko", method = "presence")
#'
NULL


#' @importFrom Matrix crossprod
#' @importFrom MultiFactor weave
#' @noRd
S7::method(mapModules, MultiFactor::MultiFactor) <- function(
    map, x, x.type = "uniref90", x.stratified = FALSE,
    method = c("sum", "count", "coverage", "presence"),
    verbose = TRUE, coverage.threshold = 0.8){

    # Check arguments
    if( !is.logical(verbose) ){
        stop("'verbose' should be TRUE or FALSE.", call. = FALSE)
    }
    if(coverage.threshold <= 0L | coverage.threshold > 1L ) {
        stop("'coverage.threshold' should be between 0-1.", call. = FALSE)
    }
    if(!x.type %in% colnames(map)) {
        stop("'x.type' must be found in 'map' argument. 'map' contained:\n",
        {paste(colnames(map), collapse = ", ")}, ".")
    }
    if(!isFALSE(x.stratified) && !is.character(x.stratified)) {
        stop("'x.stratified' should a character or FALSE.")
    }

    # Organise parameters
    method <- match.arg(arg = method, c("sum", "count", "coverage", "presence"))

    m2c_ind <- grep("_component$", names(map))
    c2c_ind <- grep("^component2.*_complex$", names(map))
    c2f_ind <- grep("_complex2", names(map))

    module_name  <- names(map[[m2c_ind]])[[1L]]
    complex_name <- names(map[[c2f_ind]])[[1L]]
    feature_name <- names(map[[c2f_ind]])[[2L]]
    all_modules  <- levels(map[[m2c_ind]][[1L]])


    complex2xtype <- if(feature_name == x.type) { map[[c2f_ind]] } else {
        weave(map, .by = c(complex_name, x.type))
        }

    mod2x <- weave(
        map, .by = c(module_name, feature_name)
    )
    c2c         <- as.matrix(map[[c2c_ind]], terms = c(2L, 1L))
    map_m2c <- map[[m2c_ind]]

    # Check for presence of all complex components
    if(isFALSE(x.stratified)) {
        complex2x <-  .map_AND_complex_unstrat(x, complex2xtype)

            component2x <- Matrix::crossprod(c2c, complex2x, boolArith = TRUE)
            # Assess module coverage
            module_coverage <- .map_coverage(map_m2c, component2x)

            if( verbose ){
                message(NROW(x), " features were mapped to ",
                        NROW(module_coverage), " modules.")
            }

            if( method == "coverage" ) {
                rownames(module_coverage) <- levels(map_m2c[[1L]])
                return(module_coverage)
            }
            module_present <- module_coverage >= coverage.threshold
            if( method == "presence" ) {
                rownames(module_present) <- levels(map_m2c[[1L]])
                return(module_present)
            }

            mod2x <- weave(map, .by = c(x.type, module_name), out.format = "matrix")

            if( method == "sum" ) {
                sum_table <- Matrix::crossprod(mod2x, x) * module_present

                rownames(sum_table) <- levels(map_m2c[[1L]])
                return(sum_table)
            }

            if( method == "count" ) {
                count_table <- Matrix::crossprod(mod2x, x != 0L) * module_present

                rownames(count_table) <- levels(map_m2c[[1L]])
                return(count_table)
            }
        }


    if(is.character(x.stratified)) {

        x.str  <- .x_strat(row.names(x), x.stratified, x.type)
        x.keep <- x.str[[1L]] %in% levels(complex2xtype[[x.type]])
        x <- x[ x.keep, ]
        x.ids     <- factor(x.str[[1L]][x.keep])
        row.names(x) <- x.ids
        x.subtype <- factor(x.str[[2L]][x.keep])


        complex2xtype[[x.type]] <- factor(
            complex2xtype[[x.type]], levels = levels(x.ids)
        )
        complex2xtype <- complex2xtype[!is.na(complex2xtype[[x.type]]), ]

        map_mat <- as.matrix(
            complex2xtype, terms = c(2L, 1L)
            )
        complex2x <-  lapply(split.data.frame(x, x.subtype),
                                FUN = .map_AND_complex_strat,
                                  map_mat = map_mat, x.ids = x.ids)

        component2x <- lapply(
            X = complex2x,
            function(x) Matrix::crossprod(c2c, x, boolArith = TRUE)
            )

    # Assess module coverage
    module_coverage <- lapply(component2x, .map_coverage, x = map_m2c)

    if( method == "coverage" ) {

        out <- do.call(rbind, module_coverage)
        row.names(out) <- paste(
            rep(all_modules, times = length(names(module_coverage))),
            rep(names(module_coverage), each = length(all_modules)),
            sep = x.stratified)
        return(out)
    }
    module_present <- lapply(module_coverage,  `>=`,  coverage.threshold)
    if( method == "presence" ) {
        out <- do.call(rbind, module_present)
        row.names(out) <- paste(
            rep(all_modules, times = length(names(module_coverage))),
            rep(names(module_coverage), each = length(all_modules)),
            sep = x.stratified)
        return(out)
    }

    mod2x <- weave(map, .by = c(x.type, module_name), out.format = "matrix")

    if( method %in% c("sum", "count" )) {
        res_table <- lapply(
            split.data.frame(x, x.subtype),
            .crossprod_strat, mod2x, x.ids, method
            )
        res_table <- mapply(`*`, res_table, module_present)

        out <- do.call(rbind, res_table)
        row.names(out) <- paste(
            rep(all_modules, times = length(names(module_coverage))),
            rep(names(module_coverage), each = length(all_modules)),
            sep = x.stratified)
        return(out)

    }

}
}

# map = module limkmap, with complex on the left (1L), features on the right (2L).
# x = feature table
# return = bool matrix indicating full complex presence, row, per sample, col.
#
.map_AND_complex_unstrat <- function(x, map) {
    map_mat <- as.matrix(map, terms = c(2L, 1L))
    Matrix::Matrix(
        Matrix::crossprod(map_mat, x !=0L ) >= Matrix::colSums(map_mat),
        sparse = TRUE
    )
}

# map_mat = module linkmap, with complex on the left (1L), features on the right (2L).
# x = feature table
# x.ids = x.ids
# return = bool matrix indicating full complex presence, row, per sample, col.
#
.map_AND_complex_strat <- function(x, map_mat, x.ids) {
    map_mat <- map_mat[match(row.names(x), levels(x.ids)),]
    if( prod(dim(map_mat)) < 2L ) { return(NULL) }

    Matrix::Matrix(
        Matrix::crossprod(map_mat, x !=0L ) >= Matrix::colSums(map_mat),
        sparse = TRUE, dimnames = list(colnames(map_mat), colnames(x))
    )
}

.crossprod_strat <- function(x, mod2x, x.ids, method) {
    if(method == "count") x <- x != 0L
    mod2x <- mod2x[match(row.names(x), levels(x.ids)),]
    if( prod(dim(mod2x) ) < 2L ) { return(NULL) }

    Matrix::crossprod(mod2x, x)
}

# x = module linkmap, with module on the left (1L), component on the right (2L).
# y = feature table
# return = numeric matrix indicating module coverage as prop, row, per sample, col.
#
.map_coverage <- function(x, y) {
    x_mat <- as.matrix(x, terms = c(2L, 1L))
    Matrix::Matrix(
        Matrix::crossprod(x_mat, y!=0L) / Matrix::colSums(x_mat),
        sparse = TRUE
    )
}

# x = stratified feature table, with row.names format feature<sep>subtype
# sep = separator character.
# cn colnames of returned data.frmae
# returns a LinkMap-formatted data.frame.
.x_strat <- function(x, sep, feature_name = "feature") {
    scan(
        text = x, what = list(
            feature_name = character(),
            subtype = character()
            ),
        sep = sep, quiet = TRUE
    )
}




