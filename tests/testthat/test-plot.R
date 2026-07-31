
test_that("plot", {
    
    expect_error(plotPath(ec ~ ko))
    
    expect_error(
        plotPath(graph, ec ~ ko, focus = "wrong"),
        "'focus' must be TRUE or FALSE."
    )
    
    expect_error(
        plotPath(graph, prune = TRUE),
        "'prune' must be FALSE when 'by' is not defined."
    )
    
    p <- plotPath(graph)
    pdata <- ggplot2::ggplot_build(p)$data

    expect_in(pdata[[1]]$edge_colour, "grey80")
    expect_in(pdata[[1]]$edge_alpha, TRUE)
    expect_identical(nrow(pdata[[2]]), length(graph))

    p <- plotPath(graph, ec ~ ko, prune = TRUE)
    pdata <- ggplot2::ggplot_build(p)$data
    
    expect_contains(pdata[[1]]$edge_colour, c("grey80", "red"))
    expect_contains(pdata[[1]]$edge_alpha, c(FALSE, TRUE))
    expect_identical(nrow(pdata[[2]]), 2L)
    
    nodes <- sum(!is.na(igraph::V(graph)$KEGG) | !is.na(igraph::V(graph)$Rhea))
    
    p <- plotPath(graph, res.name = c("KEGG", "Rhea"))
    pdata <- ggplot2::ggplot_build(p)$data
    
    expect_in(pdata[[1]]$edge_colour, "grey80")
    expect_contains(pdata[[1]]$edge_alpha, c(FALSE, TRUE))
    expect_identical(nrow(pdata[[2]]), nodes)
    
    p <- plotPath(graph, res.name = c("KEGG", "Rhea"), focus = TRUE)
    pdata <- ggplot2::ggplot_build(p)$data
    
    expect_in(pdata[[1]]$edge_colour, "grey80")
    expect_in(pdata[[1]]$edge_alpha, TRUE)
    expect_identical(nrow(pdata[[2]]), nodes)
})
