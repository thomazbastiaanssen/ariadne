
test_that("sparql", {
    
    expect_error(
        .get_batches(seq(100), 20, 2, 2),
        "Query limit was reached (25 > 20). Increase 'factor', 'batch.size' or 'workers' and try again."
    )
    
    ranges <- .get_batches(seq(100), 20, 2, 6)
    
    expect_length(ranges, 5L)
    
    out <- .querySPARQL(
        "uniprotkb", "GeneID", "UniProt", c("P0DTC2", "P60507", "Q9Y261"), 1e6
    )
    
    ids <- as.numeric(sub(".*/", "", out$GeneID))
    
    expect_in(ids, c(43740568, 105373297, 3170))
})