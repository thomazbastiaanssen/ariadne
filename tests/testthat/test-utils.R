test_that("utils", {
    
    values <- list(a = c(1, 3), b = c(1, 2, 4), c = 2)
    
    expect_no_error(linkmap <- as.linkmap(values))
    
    expect_contains(unique(linkmap$y), seq(4))
    expect_equal(dim(linkmap), c(6, 2))
    
    keys <- names(values)
    values <- unname(values)
    
    expect_no_error(
        linkmap <- as.linkmap(values, keys, col.names = c("keys", "values"))
    )
    
    expect_equal(names(linkmap), c("keys", "values"))
    expect_equal(unique(linkmap$keys), keys)
})