
ModuleDatabases <- list(
    GMM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GMMs.v1.07.txt",
    GBM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GBMs.v1.0.txt"
)

MappingDatabases <- list(
    ChocoPhlAn = list(
        repo = "https://zenodo.org/records/17100034/files/",
        from = c("eggnog", "go", "ko", "level4ec"),
        to = c("uniref50", "uniref90"),
        path = function(repo, from, to){
            paste0(repo, "map_", from, "_", to, ".txt.gz")
        }
    ),
    WoL = list(
        repo = "https://ftp.microbio.me/pub/wol-20April2021/",
        from = c("eggnog", "go", "ko", "orthodb", "refseq"),
        to = "uniref90",
        path = function(repo, from, to){
            paste0(repo, "function/", from, "/", from, ".map.xz")
        }
    )
)

.get_repo <- function(x){
    
    db <- MappingDatabases[[x]]
    
    combos <- expand.grid(db$from, db$to)
    rownames(combos) <- interaction(combos, sep = "2")
    
    map <- apply(combos, 1, function(row) {
        df <- data.frame(factor(), factor())
        names(df) <- row
        return(df)
    })
    
    mf <- MultiFactor(lapply(map, LinkMapDB))
    return(mf)
}
