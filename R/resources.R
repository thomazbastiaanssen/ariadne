
#' @name resources
#' @rdname resources

meta <- list(
    ChocoPhlAn = "https://zenodo.org/records/17100034/files/",
    GM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/",
    GO = "https://current.geneontology.org/ontology/external2go/",
    KEGG = "https://www.genome.jp/kegg/",
    TIGRFAMs = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/",
    UniProt = "https://sparql.uniprot.org/",
    WoL = "https://ftp.microbio.me/pub/wol-20April2021/"
)

meta <- data.frame(
  name = names(meta),
  repo = unlist(meta, use.names = FALSE)
)
