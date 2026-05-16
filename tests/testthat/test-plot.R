
test_that("plot", {
    
    graph <- ariadne()
    
    expect_error(plotPath(ec ~ ko))
    
    p <- plotPath(graph)
    pdata <- ggplot2::ggplot_build(p)$data

    expect_in(pdata[[1]]$edge_colour, "grey80")
    expect_in(pdata[[1]]$edge_alpha, TRUE)
    expect_identical(nrow(pdata[[2]]), length(graph))

    p <- plotPath(graph, ec ~ ko, prune = TRUE)
    pdata <- ggplot2::ggplot_build(p)$data
    
    expect_in(pdata[[1]]$edge_colour, c("grey80", "red"))
    expect_in(pdata[[1]]$edge_alpha, c(FALSE, TRUE))
    expect_identical(nrow(pdata[[2]]), 2L)
    
    expect_error(plotPath(graph, ec ~ ko, focus = "wrong"))
    
    expect_error(
        plotPath(graph, prune = TRUE),
        "'prune' must be FALSE when 'by' is not defined."
    )
})
