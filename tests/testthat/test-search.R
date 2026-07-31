test_that("search", {
    
    expect_error(searchPath(ko ~ ec))
  
    expect_message(
        expect_invisible(
            searchPath(graph, ko ~ ec)
        )
    )
})