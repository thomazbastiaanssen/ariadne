test_that("importMapping", {
  
    expect_error(
        importMapping("ChocoPhlAn", from = "wrong", to = "uniref90"),
        "'from' should be defined and be one of eggnog, go, ko, level4ec."
    )
    
    expect_error(
        importMapping("ChocoPhlAn", from = "ko", to = "wrong"),
        "'to' should be defined and be one of uniref50, uniref90."
    )
    
    expect_no_error(
        map <- importMapping("ChocoPhlAn", from = "ko", to = "uniref90")
    )
    
    expect_named(map)
  
    expect_error(
        importMapping(
            "ChocoPhlAn",
            from = c("eggnog", "ko"),
            to = c("uniref90", "uniref90", "uniref90")
        ),
        "'map.files', 'from' and 'to' must have compatible lengths."
    )
    
    expect_no_error(
        importMapping("ChocoPhlAn", from = c("eggnog", "ko"), to = "uniref90")
    )
})
