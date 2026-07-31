
test_that("custom", {
    
    url <- "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv"
    
    expect_error(addResource(graph, url, select = c("RHEA_ID", "ID"),
        col.names = c("x", "y")),
        "At least one feature must be in 'graph'.", fixed = TRUE
    )
    
    expect_error(addResource(graph, url, res.name = "Rhea",
        select = c("RHEA_ID", "ID"), col.names = c("rhea", "ec")),
        "This link for Rhea already exists. Set 'force' to TRUE to overwrite it.",
        fixed = TRUE
    )
    
    edge_df <- graph |>
        addResource(
            url, res.name = "Rhea", force = TRUE,
            select = c("RHEA_ID", "ID"), col.names = c("rhea", "ec")
        ) |>
        igraph::as_data_frame(what = "edges")
    
    expect_contains(edge_df$url, url)
    expect_equal(sum(.get_edge_keys(edge_df) == "ec_rhea_Rhea"), 1L)
})