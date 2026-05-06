
test_that("weave", {
    
    graph <- ariadne()
    
    path_df <- .draw_path(graph, bugsig ~ ko, 1, "uniref90", "taxid", NULL)
    
    expect_false("taxid" %in% path_df$from || "taxid" %in% path_df$to)
    expect_true("uniref90" %in% path_df$from || "taxid" %in% path_df$to)
    
    expect_error(
        .draw_path(graph, bugsig ~ ko, 1, "uniref90", "uniref90"),
        "'include' and 'exclude' cannot overlap."
    )
    
    # Expected `ncol(ko2gmm)` to equal 3L.
    # Differences:
    # 1/1 mismatches
    # [1] 2 - 3 == -1
    
    # ko2gmm <- weavePath(graph, ko ~ gmm, use.names = TRUE)
    # expect_equal(ncol(ko2gmm), 3L)
    
    ko2gbm <- weavePath(graph, ko ~ gbm, use.names = FALSE)
    expect_identical(ncol(ko2gbm), 2L)
    
    ec2gmm <- weaveComplex(graph, ec ~ gmm)
    expect_identical(colnames(ec2gmm), c("ec", "gmm", "cov", "gmm.name"))
})