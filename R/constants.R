
ChocoPhlAn <- list(
    eggnog.uniref50 = "map_eggnog_uniref50.txt.gz",
    eggnog.uniref90 = "map_eggnog_uniref90.txt.gz",
    go.uniref50 = "map_go_uniref50.txt.gz",
    go.uniref90 = "map_go_uniref90.txt.gz",
    ko.uniref50 = "map_ko_uniref50.txt.gz",
    ko.uniref90 = "map_ko_uniref90.txt.gz",
    ec.uniref50 = "map_level4ec_uniref50.txt.gz",
    ec.uniref90 = "map_level4ec_uniref90.txt.gz"
)

ChocoPhlAn <- lapply(ChocoPhlAn, function(x)
    paste0("https://zenodo.org/records/17100034/files/", x))

GutModules <- list(
    GMM = "GMMs.v1.07.txt",
    GBM = "GBMs.v1.0.txt"
)

GutModules <- lapply(GutModules, function(x)
    paste0("https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/", x))
