test_that("getModules", {
    
    data("butyrate", package = "ariadne")
    tax.before <- nrow(butyrate)
    
    butyrate <- getFullTaxonomyLabels(butyrate)
    
    expect_false(is.data.frame(butyrate))
    tax.after <- length(unlist(butyrate))
    
    expect_equal(tax.before, tax.after)
    
    data("Tengeler2020", package = "mia")
    tse <- Tengeler2020
    
    modules <- getModules(tse, butyrate)
    tse <- addModules(tse, butyrate)
    
    expect_length(modules, nrow(tse))
    
    expect_error(
        getModules(tse, butyrate, by = 2L),
        "'group' does not match any variables in the side information."
    )
    
    modules <- list(A = c("g__Bacteroides", "f__Enterobacteriaceae"),
                    B = c("c__Clostridia", "p__Cyanobacteria"))
    
    expect_error(
        getModules(tse, modules, exact.tax.level = "wrong"),
        "'modules' must be a character vector or list of character vectors, where each vector corresponds to a module."
    )
    
    mod.table <- getModules(tse, modules, exact.tax.level = TRUE)
    
    col.modules <- list(A = c("Cohort_1", "Cohort_3"),
                        B = c("Cohort_2"),
                        c("Cohort_1", "Cohort_2", "Cohort_3"))
    
    expect_error(
        getModules(tse, col.modules, by = 2L),
        "'group' does not match any variables in the side information."
    )
    
    mod.table <- getModules(tse, col.modules, by = 2L, group = "cohort")
})