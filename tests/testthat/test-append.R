
test_that("append", {
    
    data("Tengeler2020", package = "mia")
    data("butyrate", package = "ariadne")
    
    tse <- Tengeler2020
    
    expect_error(
        getModules(tse, butyrate, key = "wrong"),
        "'key' must be 'row.names' or a variable of 'x'.",
        fixed = TRUE
    )
    
    expect_error(
        getModules(tse, butyrate, by = "wrong"),
        "'by' must be either 'rows' or 'cols'.",
        fixed = TRUE
    )
    
    expect_error(
        getModules(tse, butyrate, key = "Genus", as = "wrong"),
        "'as' must be either 'ids' or 'names'",
        fixed = TRUE
    )
    
    col.names <- c(names(rowData(tse)), "butyrate")
    
    modules <- getModules(tse, butyrate, key = "Genus")
    tse <- addModules(tse, butyrate, key = "Genus", as = "names")
    
    expect_s3_class(modules, "data.frame")
    expect_in(colnames(modules), "C00246")
    expect_in(rownames(modules), rownames(tse))
    expect_contains(modules$C00246, c(TRUE, FALSE))
    
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_in(names(rowData(tse)), col.names)
    expect_contains(rowData(tse)$butyrate, c(TRUE, FALSE))
})