test_that("humann", {
    
    ko_ids <- c(
        # MF0002
        "K03332",
        # MF0003
        "K01051", "K01184", "K01213",
        # MF0004
        "K01731", "K01625",
        # MF0007
        "K02786", "K01635"
    )
    
    # ko2uniref <- weavePath(graph, ko ~ uniref90, init = ko_ids)
    # ko2uniref <- ko2uniref[!duplicated(ko2uniref$ko), ]
    
    uniref_ids <- c(
        "UniRef90_A0A010Q3S3", "UniRef90_A0A010QND4", "UniRef90_A0A010QRC3",
        "UniRef90_A0A023WSR0", "UniRef90_A0A060B5N9", "UniRef90_A0A068TIR7",
        "UniRef90_A0A075LLW4", "UniRef90_A0A0E1NTX0"
    )
    
    uniref2ko <- weavePath(graph, uniref90 ~ ko, init = uniref_ids)
    
    tax2uniref <- data.frame(taxname = "ariadne bug", uniref90 = uniref_ids)
    
    tax2gmm <- weaveComplex(graph, taxname ~ gmm, init = tax2uniref)
    
    n_features <- length(uniref_ids)
    n_samples <- 3
    max_abund <- 5
    
    set.seed(123)
    
    mat <- seq(max_abund) |>
        sample(n_features * n_samples, replace = TRUE) |>
        matrix(n_features, n_samples)
    
    rownames(mat) <- uniref_ids
    colnames(mat) <- letters[1:n_samples]
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat)
    )
    
    rowData(se)$taxname <- "ariadne bug"
    
    ko_mat <- se |>
        addModules(uniref2ko, "rows") |>
        mia::agglomerateByModule(1L, uniref2ko$ko) |>
        SummarizedExperiment::assay()
    
    gmm_mat <- se |>
        mia::agglomerateByVariable(1L, "taxname") |>
        addModules(tax2gmm, "rows") |>
        mia::agglomerateByModule(1L, tax2gmm$gmm) |>
        SummarizedExperiment::assay()
    
    #expect_identical()
    #expect_equal()
    # against humann and omixer outputs

})