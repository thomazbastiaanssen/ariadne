
test_that("plot", {
    
    graph <- ariadne()
    
    expect_error(plotPath(ec ~ ko))
    
    expect_no_error(plotPath(graph, ec ~ ko))
    
    expect_error(plotPath(graph, ec ~ ko, focus = "wrong"))

})