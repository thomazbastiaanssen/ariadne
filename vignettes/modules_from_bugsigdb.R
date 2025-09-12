# Import libraries
devtools::load_all() # On mia/module branch
library(bugsigdbr)

# Download database
bsdb <- importBugSigDB()

# Retrieve microbial signatures
sigs <- getSignatures(bsdb, tax.id.type = "metaphlan")

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

# Add modules table to rowData
tse <- addModules(tse, sigs)
# Agglomerate features by ADHD modules
tse.mod <- agglomerateByModule(tse, by = 1L, group = names(sigs))

# Select non-empty modules
keep <- rowSums(assay(tse.mod)) > 0
tse.mod <- tse.mod[keep, ]



data("butyrate", package = "ariadne")

butyrate

sigs <- .tax_table2label(butyrate)

tse <- addModules(tse, sigs)

