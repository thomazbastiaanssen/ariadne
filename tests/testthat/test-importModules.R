test_that("importModules", {

    expect_warning(
        expect_error(
            importModules("wrong")
        )
    )

    expect_no_error(gbm <- importModules("GBM"))

})
