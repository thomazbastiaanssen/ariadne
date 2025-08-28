
# 3 options (module.file, table, list)
# addModules(tse, sigs = "GBM", ids = rowdata.colname)
# addModules(tse, sigs = butyrate)
# addModules(tse, sigs = sig.list)

gbm <- importModules("GBM")

map <- importMapping("data/map_ko_uniref90.txt")

sigs <- mapModules(gbm, map)

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

modules <- getModules(tse, sigs)

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))



# The following approach works

qls <- GET("https://rest.uniprot.org/uniref/UniRef90_A0A061EXN5/members?format=list")
content_text <- content(qls, as = "text")
uniref_ids <- unlist(strsplit(content_text, "\n"))

first_entry <- queryUniProt(
  query = uniref_ids,
  fields = c("id", "organism_name", "lineage"),
  collapse = " OR "
)
