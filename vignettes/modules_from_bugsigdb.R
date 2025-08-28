# Import libraries
devtools::load_all() # On mia/module branch
library(bugsigdbr)

# Download database
bsdb <- importBugSigDB()

names(bsdb)

# Select studies related to ADHD
adhd.bugs <- bsdb |>
    subset(grepl("ADHD", Title))

# Retrieve signatures from selected studies
adhd.sigs <- getSignatures(
    adhd.bugs,
    tax.id.type = "metaphlan"
)

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

# Build modules adjacency table from bugsigdb signatures
modules <- getModules(tse, adhd.sigs)
# Add modules table to rowData
tse <- addModules(tse, adhd.sigs, exact.tax.level = FALSE)
# Agglomerate features by ADHD modules
tse.mod <- agglomerateByModule(tse, by = 1L, group = names(adhd.sigs))

# save(butyrate, file = "data/butyrate.rda")
# butyrate <- read.csv("data/butyrate.csv")

data("butyrate", package = "mia")

butyrate

sigs <- .tax_table2label(butyrate)

tse <- addModules(tse, sigs)

library(biomaRt)

data("GlobalPatterns")
tse <- GlobalPatterns

mart <- useMart("ENSEMBL_MART_ENSEMBL")

dataset <- useDataset("uniprotswissprot", mart)

protein_names <- getBM(attributes = c("uniprotswissprot"), mart = ensembl)


up <- UniProt.ws(taxId = 9606)  # Use the appropriate tax ID for your organism
results <- select(up, keys = "UniRef90_A0A0C7KDG2",
                  columns = c("Organism", "UniRef90"), keytype = "UniRef90")

library(httr)

# Example of a UniRef90 ID
uniref_id <- "UniRef90_A4FHA7"  # Replace with your actual UniRef90 ID

response <- GET(paste0("https://www.uniprot.org/uniref/", uniref_id))
content <- content(response, "text")

library(rvest)

GET("https://rest.uniprot.org/uniref/UniRef90_A0A061EXN5/members?format=list&size=500")

uniref_id <- "UniRef90_A0A061EXN5"
link <- paste0("https://www.uniprot.org/uniref/", uniref_id)

page <- read_html(link)

page |>
  html_elements(".S2LRp a") |>
  html_text()
  
mapUniProt(from = "UniRef90", to = "Taxonomy",
           query = "UniRef90_A4FHA7")

uniref_ids <- c("UniRef90_A0A059A5P0", "UniRef90_A0A061EYG9", "UniRef90_A0A068NE04",
  "UniRef90_A0A068T9S4", "UniRef90_A0A077LNQ7")

ko_uniref90_dict <- read.table("~/Downloads/KO_uniref90_dict.txt", sep = "\t")

uniref_ids <- gsub("UniRef90_", "", uniref_ids)

parse_file <- function(file_path) {
  # Read the file content
  lines <- readLines(file_path)
  
  vectors <- list()  # Initialize an empty list to store the vectors
  names <- list()
  
  for (line in lines) {
    items <- unlist(strsplit(line, "\t"))  # Split by tab character
    vectors <- append(vectors, list(items[-1]))  # Append the entire line as a vector
    names <- append(names, items[1])
  }
  
  names(vectors) <- names
  
  return(vectors)
}

# Example usage
parsed_vectors <- parse_file("~/Downloads/KO_uniref90_dict.txt")

parsed_vectors <- lapply(parsed_vectors, function(x) gsub("UniRef90_", "", x))

library(UniProt.ws)

mapUniProt(
  from = "UniProtKB_AC-ID",
  to = "UniProtKB",
  query = parsed_vectors[[1]],
  columns = "organism_name"
)

library(KEGGREST)

query <- keggGet("K22014", "kgml")
query$ENTRY

library(xml2)
library(rvest)

link <- "https://www.kegg.jp/kegg-bin/show_brite?htext=br08611&hier=20&highlight=join_brite&mapper=taxmap%2ds%20K22014&pruning=taxmap&selected=none&extend=H1-11082"
page <- read_html(link)

text <- page |> html_elements("script")

text <- text[7]

real_text <- text |> html_text()

grep("Escherichia coli K-12 MG1655", real_text)


