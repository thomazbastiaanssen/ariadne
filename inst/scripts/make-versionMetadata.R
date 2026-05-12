
devtools::load_all()
source("utils.R")

default.graph <- 19397292
dfs <- list()

# Add BugSigDB versions
record_id <- 5606166
versions <- fetch_zenodo_versions(record_id)

dfs[["BugSigDB"]] <- data.frame(
    version = names(versions), key = versions, graph = default.graph
)

# Add ChocoPhlAn versions
record_id <- 17100034
versions <- fetch_zenodo_versions(record_id)

dfs[["ChocoPhlAn"]] <- data.frame(
    version = names(versions), key = versions, graph = default.graph
)

# Add GO versions
versions <- c("2026-03-25", "2026-01-23", "2025-10-10")

dfs[["GO"]] <- data.frame(
    version = versions, key = versions, graph = default.graph
)

# Add WoL versions
dfs[["WoL"]] <- data.frame(
    version = c("v2", "v20April2021"),
    key = c("wol2", "wol-20April2021"),
    graph = c(default.graph, 18788726)
)

# Add Misc versions
dfs[["GM"]] <- data.frame(
    version = "v1",
    key = "omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata",
    graph = default.graph
)

# Add TIGRFAMs versions
dfs[["TIGRFAMs"]] <- data.frame(
    version = "v15", key = "release_15.0", graph = default.graph
)

# Add unversioned resources
not_versioned <- c("KEGG", "OTT", "Rhea", "UniProt")
nover_dfs <- lapply(
    not_versioned,
    function(x) data.frame(version = "latest", key = "", graph = default.graph)
)
# Add names to dataframes for unversioned resources
names(nover_dfs) <- not_versioned

# Bind all resource dataframes
df <- bind_rows(c(dfs, nover_dfs), .id = "source")
rownames(df) <- NULL

# Set latest versions as default
df$default <- !duplicated(df$source)
# Order resources alphabetically
df <- df[order(df$source), ]

# default.versions <- c(
#     BugSigDB = "v1.3.0", ChocoPhlAn = "v201901b", Misc = "v1",
#     GO = "2026-03-25", KEGG = NA, OTT = NA, Rhea = NA, TIGRFAMs = "v15",
#     UniProt = NA, WoL = "v2"
# )

# Add resource base urls as attribute
attr(df, "urls") <- c(
    BugSigDB = "https://zenodo.org/records/",
    ChocoPhlAn = "https://zenodo.org/records/",
    GM = "https://github.com/",
    GO = "https://release.geneontology.org/",
    KEGG = "https://www.genome.jp/kegg",
    OTT = "https://opentreeoflife.github.io",
    Rhea = "https://www.rhea-db.org",
    TIGRFAMs = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/",
    UniProt = "https://www.uniprot.org",
    WoL = "https://ftp.microbio.me/pub/"
)

df -> versionMetadata
usethis::use_data(versionMetadata, internal = TRUE, overwrite = TRUE)
