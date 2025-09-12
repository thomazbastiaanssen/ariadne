# 3 options (module.file, table, list)
# addModules(tse, sigs = "GBM", ids = rowdata.colname)
# addModules(tse, sigs = butyrate)
# addModules(tse, sigs = sig.list)



mod.ko <- importModules("GMM")
ko.uniref90 <- importMapping("ko.uniref90")

sigs <- mapModules(mod.ko, ko.uniref90, drop.unmatched = TRUE)

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

modules <- getModules(tse, sigs)

tse <- addModules(tse, sigs)

tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))
