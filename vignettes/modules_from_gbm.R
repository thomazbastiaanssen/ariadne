
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

sigs <- Filter(function(sig) length(sig) > 0, sigs)

modules <- getModules(tse, sigs)

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))
