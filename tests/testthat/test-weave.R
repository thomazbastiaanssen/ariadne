
test_that("weave", {
    
    graph <- ariadne()
    
    path_df <- .draw_path(graph, bugsig ~ ko, 1, "uniref90", "taxid")
    
    expect_false("taxid" %in% path_df$from || "taxid" %in% path_df$to)
    expect_true("uniref90" %in% path_df$from || "taxid" %in% path_df$to)
    
    expect_error(
        .draw_path(graph, bugsig ~ ko, 1, "uniref90", "uniref90"),
        "'include' and 'exclude' cannot overlap."
    )
    
    expect_warning(
        weavePath(graph, ko ~ ec, k = 1, use.names = TRUE),
        "Names for 3 ec ids not found."
    )
    
    
    #weaveComplex(graph,  ~ gmm)
    
})