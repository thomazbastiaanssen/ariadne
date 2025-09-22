test_that("mapModules", {
    
    gbm <- importModules("GBM")
    eggnog.uniref90 <- importMapping("ChocoPhlAn", "eggnog", "uniref90")
    
    expect_error(
        mapModules(eggnog.uniref90, gbm),
        "'map' did not match any element in 'modules'."
    )
    
    expect_error(
        mapModules(gbm, eggnog.uniref90, mode = "wrong"),
        "'mode' must be either uniref or taxonomy."
    )
    
    sig.list1 <- mapModules(
        head(gbm),
        eggnog.uniref90,
        mode = "taxonomy",
        type = "andor",
        remove.empty = TRUE)
    
    sig.list2 <- mapModules(
        head(gbm),
        eggnog.uniref90,
        mode = "taxonomy",
        type = "andor",
        remove.empty = FALSE
    )
    
    expect_length(sig.list1, 2)
    expect_length(sig.list2, 6)
    
})
