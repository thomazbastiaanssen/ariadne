# 3 options (module.file, table, list)
# addModules(tse, sigs = "GBM", ids = rowdata.colname)
# addModules(tse, sigs = butyrate)
# addModules(tse, sigs = sig.list)

gbm <- importModules("GBM")

map1 <- importMapping("ko.uniref90")
map2 <- importMapping("eggnog.uniref90")
map3 <- importMappings(c("ko.uniref90", "eggnog.uniref90"))

sigs <- mapModules(gbm[seq(3)], map1)

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

modules <- getModules(tse, sigs)

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))
