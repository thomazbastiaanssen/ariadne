
# 3 options (module.file, table, list)
# addModules(tse, sigs = "GBM", ids = rowdata.colname)
# addModules(tse, sigs = butyrate)
# addModules(tse, sigs = sig.list)

gbm <- importModules("GBM")

map <- importMapping("data/map_ko_uniref90.txt")

sigs <- mapModules(gbm, map) # reasonably fast


# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

modules <- getModules(tse, sigs) # quite slow

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))
