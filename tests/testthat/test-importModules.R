test_that("importModules", {
    
    expect_error(importModules("wrong"))
    
    expect_no_error(gbm <- importModules("GBM"))
    expect_no_error(gms <- importModules(c("GBM", "GMM")))
    
    expect_length(gbm, 56)
    expect_length(gms, 159)
})