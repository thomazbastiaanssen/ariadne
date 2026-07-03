test_that("utils", {
    
    expect_error(
        listResourceVersions(default = c(NA, NA)),
        "'default' must be TRUE or FALSE."
    )
    
    meta <- listResourceVersions()
    default <- listResourceVersions(default = TRUE)
    
    expect_gt(anyDuplicated(meta$resource), 0L)
    expect_equal(anyDuplicated(default$resource), 0L)
    
    
    set.seed(123)
    
    gene_ids <- c(
        "Unbinned",
        "uniref90|Bacteroides.thetaiotaomicron",
        "uniref90|Escherichia.coli",
        "uniref90|Faecalibacterium.prausnitzii",
        "uniref90|unclassified",
        "uniref90|Bifidobacterium.adolescentis",
        "uniref90|Lactobacillus.rhamnosus"
    )
    
    # Create a matrix of random counts (e.g., 6 genes x 10 samples)
    count_matrix <- matrix(
        data = sample(10:1000, size = 10 * length(gene_ids), replace = TRUE),
        nrow = length(gene_ids)
    )
    
    rownames(count_matrix) <- gene_ids
    
    # Create the SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = count_matrix)
    )
    
    se <- processGeneFamilies(se)
    
    expect_s4_class(se, "SummarizedExperiment")
    expect_contains(
        names(SummarizedExperiment::rowData(se)),
        c("uniref90", "taxname", "genus", "species")
    )
})