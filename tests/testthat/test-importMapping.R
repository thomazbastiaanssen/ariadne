test_that("importMapping", {

    expect_error(
        importMapping("ChocoPhlAn", subset = c("wrong", "uniref90")),
        "both terms must be found as colnames in 'x'"
    )

  expect_error(
    importMapping("ChocoPhlAn", subset = c("ko", "wrong")),
    "both terms must be found as colnames in 'x'"
  )

    map <- importMapping("ChocoPhlAn", subset = c("ko", "uniref90"),
                         dry_run = FALSE)

    expect_named(map)


})
