
#' @name resources
#' @rdname resources

.create_resource <- function(repo = NULL, version = NULL, from = NULL,
    to = NULL, path = NULL){
    # Add dummy names to unnamed features
    if(is.null(names(from))){
        names(from) <- rep("", length(from))
    }
    if(is.null(names(to))){
        names(to) <- rep("", length(to))
    }
    # Add proper names to features
    names(from) <- ifelse(names(from) == "", unlist(from), names(from))
    names(to) <- ifelse(names(to) == "", unlist(to), names(to))
    # Create resource
    resource <- list(
        repo = repo,
        version = version,
        from = from,
        to = to,
        path = path
    )
    return(resource)
}

Resources <- list(
    ChocoPhlAn = .create_resource(
        repo = "https://zenodo.org/records/",
        version = "17100034/files/",
        from = c("eggnog", "go", "ko", ec = "level4ec"),
        to = c("uniref50", "uniref90"),
        path = function(repo, version, from, to){
            paste0(repo, version, "map_", from, "_", to, ".txt.gz")
        }
    ),
    GMM = .create_resource(
        repo = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/",
        version = "GMMs.v1.07.txt",
        from = "gmm",
        to = "ko",
        path = function(repo, version, from, to){
            paste0(repo, version)
        }
    ),
    GBM = .create_resource(
        repo = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/",
        version = "GBMs.v1.0.txt",
        from = "gbm",
        to = c("eggnog", "ko", "tigr"),
        path = function(repo, version, from, to){
            paste0(repo, version)
        }
    ),
    GO = .create_resource(
        repo = "https://current.geneontology.org/ontology/external2go/",
        from = c("ec", "hamap", "interpro", reaction = "kegg_reaction",
            "metacyc", "pfam", "pirsf", "prints", "prosite", "reactome",
            "resid", "rfam", "rhea", "smart", "um-bbd_enzymeid",
            "um-bbd_pathwayid", "um-bbd_reactionid", "uniprotkb_kw",
            "uniprotkb_sl", "unirule", "wikipedia"),
        to = "go",
        path = function(repo, version, from, to){
            paste0(repo, version, from, "2", to)
        }
    ),
    KEGG = .create_resource(
        from = c("pathway", "brite", kegg.mod = "module", "vg", "ag",
                 "reaction", "rclass", "enzyme", "disease", "drug", "ec"),
        to = c("ko", "cpd")
    ),
    UniProt = .create_resource(
        from = c("uniref50", "uniref90"),
        to = "taxonomy"
    ),
    WoL = .create_resource(
        repo = "https://ftp.microbio.me/pub/",
        version = "wol-20April2021/",
        from = c("eggnog", "go", ko = "kegg", "orthodb", "refseq"),
        to = "uniref90",
        path = function(repo, version, from, to){
            # Account for exceptions
            suffix <- ifelse(from == "go", "all", from)
            suffix <- ifelse(from == "kegg", "ko", suffix)
            # Create path
            paste0(repo, version, "function/", from, "/", suffix, ".map.xz")
        }
    ),
    TIGRFAMs = .create_resource(
        repo = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/",
        version = "release_15.0/",
        from = c(tigr = "tigrfams"),
        to = "go",
        path = function(repo, version, from, to){
            paste0(repo, version, toupper(from), "_", toupper(to), "_LINK")
        }
    )
)
