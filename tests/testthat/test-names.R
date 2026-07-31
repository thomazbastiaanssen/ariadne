test_that("names", {
    
    expect_error(
        linkNames(graph, "wrong"), "'x' must be in 'graph'.", fixed = TRUE
    )
    
    expect_error(
        linkNames(graph, "ko", ids = "K0001", names = "alcohol dehydrogenase"),
        "Either 'ids' or 'names' can be specified.", fixed = TRUE
    )
    
    x <- linkNames(graph, "gmm")
    expect_named(x, c("gmm", "gmm.name"))
    
    name.vec <- c("alcohol dehydrogenase", "alditol oxidase")
    x <- linkNames(graph, "ec", names = name.vec)
    
    expect_equal(x[[2L]], name.vec)
    
    expect_warning(
        x <- linkNames(graph, "bugsig", ids = c("83/1/1", "wrong", "39/2/1")),
        "names for 1 bugsig ids not found.", fixed = TRUE         
    )
    
    expect_equal(is.na(x[[2L]]), c(FALSE, TRUE, FALSE))
})