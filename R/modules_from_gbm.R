devtools::load_all() # mia
library(UniProt.ws)
library(httr)
library(tidyr)

KO2UniRef90 <- importMapping("data/KO_uniref90_dict.txt")
KO2UniRef90 <- lapply(KO2UniRef90, function(x) gsub("UniRef90_", "", x))

gbm <- importModules("GMM")

sigs <- mapTaxonomy(gbm, from = "KO")

# Fix to group by mod col
sigs <- mia:::.tax_table2label(sigs[ , TAXONOMY_RANKS[-2]])


# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

mods <- getModules(tse, sigs)

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = "module")





# The following approach works

qls <- GET("https://rest.uniprot.org/uniref/UniRef90_A0A061EXN5/members?format=list")
content_text <- content(qls, as = "text")
uniref_ids <- unlist(strsplit(content_text, "\n"))

first_entry <- queryUniProt(
  query = uniref_ids,
  fields = c("id", "organism_name", "lineage"),
  collapse = " OR "
)
