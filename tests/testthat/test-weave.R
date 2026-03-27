
test_that("weave", {
    
    graph <- ariadne()
    
    path_df <- .draw_path(graph, bugsig ~ ko, 1, "uniref90", "taxid")
    
    expect_false("taxid" %in% path_df$from || "taxid" %in% path_df$to)
    expect_true("uniref90" %in% path_df$from || "taxid" %in% path_df$to)
    
    expect_error(
        .draw_path(graph, bugsig ~ ko, 1, "uniref90", "uniref90"),
        "'include' and 'exclude' cannot overlap."
    )
    
    ko2gmm <- weavePath(graph, ko ~ gmm, use.names = TRUE)
    expect_equal(ncol(ko2gmm), 3L)
    
    ko2gbm <- weavePath(graph, ko ~ gbm, use.names = FALSE)
    expect_equal(ncol(ko2gbm), 2L)
    
    expect_error(
        weaveComplex(graph, gmm ~ gbm),
        "Exactly one side of 'by' must be a module name."
    )
})