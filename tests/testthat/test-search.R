
test_that("search", {
    
    graph <- ariadne()
    
    expect_error(searchPath(ko ~ ec))
  
    expect_message(
        expect_invisible(
            searchPath(graph, ko ~ ec)
        )
    )
})