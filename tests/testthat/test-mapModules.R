test_that("mapModules", {
    
    gbm <- importModules("GBM")
    eggnog.uniref90 <- importMapping("ChocoPhlAn", eggnog ~ uniref90)
    
    expect_error(
        mapModules(eggnog.uniref90, gbm),
        "Items in 'x' did not match any item in 'y'."
    )
    
    expect_error(
        mapModules(gbm, eggnog.uniref90, mode = "wrong"),
        "'mode' must be either single or andor."
    )
    
    sig.list1 <- mapModules(
        head(gbm),
        eggnog.uniref90,
        mode = "andor",
        uniprot = TRUE,
        remove.empty = TRUE
    )
    
    sig.list2 <- mapModules(
        head(gbm),
        eggnog.uniref90,
        mode = "andor",
        uniprot = TRUE,
        remove.empty = FALSE
    )
    
    expect_length(sig.list1, 2)
    expect_length(sig.list2, 6)
})
