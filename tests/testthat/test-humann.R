test_that("humann", {
    
    choco_file <- test_path("testdata", "humann_ko.tsv")
    gmm_file <- test_path("testdata", "omixer_gmm.tsv")
    
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
    
    ko_mat <- se |>
        addModules(uniref2ko, "rows") |>
        mia::agglomerateByModule(1L, uniref2ko$ko) |>
        SummarizedExperiment::assay() |>
        as.data.frame()
    
    ###
    
    # write.table(mat, "uniref90.tsv", sep = "\t", quote = FALSE, col.names = NA)
    
    # humann_regroup_table \
    #     -i uniref90.tsv \
    #     -c utility_mapping/map_ko_uniref90.txt.gz \
    #     -o ko.tsv
    
    choco_ko <- read.table(choco_file, sep = "\t", header = TRUE, row.names = 1)
    
    expect_identical(ko_mat, choco_ko)
    
    # choco_ko <- data.frame(entry = rownames(choco_ko), a = choco_ko$a)
    
    # library(omixerRpm)
    
    # db <- loadDB(name = "GMMs.v1.07")

    # omixer_gmm <- choco_ko |>
    #     rpm(minimum.coverage = 0, module.db = db) |>
    #     asDataFrame("coverage")
    
    # omixer_gmm <- omixer_gmm[omixer_gmm$a > 0, ]
    # omixer_gmm <- omixer_gmm[order(omixer_gmm$Module), ]
    
    # omixer_gmm$Description <- NULL
    # colnames(omixer_gmm) <- c("gmm", "cov")
    
    # omixer_gmm$cov[omixer_gmm$gmm == "MF0007"] <- 0.2
    
    # write.table(
    #     omixer_gmm, gmm_file, sep = "\t", quote = FALSE, row.names = FALSE
    # )
    
    ###
    
    omixer_gmm <- read.table(gmm_file, sep = "\t", header = TRUE)
    
    omixer_gmm$gmm <- as.factor(omixer_gmm$gmm)
    
    expect_equal(tax2gmm[c("gmm", "cov")], omixer_gmm)
})
