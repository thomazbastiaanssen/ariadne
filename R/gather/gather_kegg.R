#Aliases
ko_ec     <- read.delim("https://rest.kegg.jp/link/ko/ec",    header = F, quote = "")
cpd_drug  <- read.delim("https://rest.kegg.jp/link/cpd/drug", header = F, quote = "")

#Bridges
cpd_ec    <-  read.delim("https://rest.kegg.jp/link/cpd/ec",  header = F, quote = "")
drug_ko   <-  read.delim("https://rest.kegg.jp/link/drug/ko", header = F, quote = "")

#ID Aliases
cpd_labs  <-  read.delim("https://rest.kegg.jp/list/cpd",     header = F, quote = "")
drug_labs <-  read.delim("https://rest.kegg.jp/list/drug",    header = F, quote = "")
ec_labs   <-  read.delim("https://rest.kegg.jp/list/ec",      header = F, quote = "")
ko_labs   <-  read.delim("https://rest.kegg.jp/list/ko",      header = F, quote = "")

colnames(ko_ec)        = c("ec",  "ko" )
colnames(cpd_drug)     = c("dr",  "cpd")

colnames(cpd_ec)       = c("ec",  "cpd")
colnames(drug_ko)      = c("ko",  "dr")

colnames(cpd_labs)     = c("cpd", "cpd_labs")
colnames(drug_labs)    = c("dr",  "dr_labs")
colnames(ec_labs)      = c("ec",  "ec_labs")
colnames(ko_labs)      = c("ko",  "ko_labs")


save(ko_ec, cpd_drug, 
     cpd_ec, drug_ko, 
     cpd_labs, drug_labs, ec_labs, ko_labs, file = "data/ariadne.RData")

#load("data/ariadne.RData")
