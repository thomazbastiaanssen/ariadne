# # Dummy data
# mod_tab <- importModules("GBM")
# x <- mod_tab$ko_complex2ko
# y <- replicate(20, sample(c(TRUE, FALSE), size = NROW(x), replace = TRUE), simplify = TRUE)


# x = module limkmap, with complex on the left (1L), features on the right (2L).
# y = feature table
# return = bool matrix indicating full complex presence, row, per sample, col.
#
.map_AND_complex <- function(x, y) {
    x_mat <- Matrix::sparseMatrix(
        i = x[[2L]], j = x[[1L]],
        dims = c(NROW(y), nlevels(x[[1L]]))
    )
Matrix::crossprod(x_mat, y!=0L) >= Matrix::colSums(x_mat)
}



# # Dummy data
# x <- mod_tab$GBM2GBM_component
# y <- replicate(20, sample(c(TRUE, FALSE), size = NROW(x), replace = TRUE), simplify = TRUE)


# x = module limkmap, with module on the left (1L), component on the right (2L).
# y = feature table
# return = numeric matrix indicating module coverage as prop, row, per sample, col.
#
.map_coverage <- function(x, y) {
    x_mat <- Matrix::sparseMatrix(
        i = x[[2L]], j = x[[1L]],
        dims = c(NROW(y), nlevels(x[[1L]]))
    )
    Matrix::crossprod(x_mat, y!=0L) / Matrix::colSums(x_mat)
}
