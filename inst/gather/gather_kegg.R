#Aliases
ko_ec       <- read.delim("https://rest.kegg.jp/link/ko/ec",    header = F, quote = "")
cpd_drug    <- read.delim("https://rest.kegg.jp/link/cpd/drug", header = F, quote = "")
cpd_glycan  <- read.delim("https://rest.kegg.jp/link/cpd/glycan", header = F, quote = "")

#Bridges
cpd_ec    <-  read.delim("https://rest.kegg.jp/link/cpd/ec",  header = F, quote = "")
drug_ko   <-  read.delim("https://rest.kegg.jp/link/drug/ko", header = F, quote = "")
glycan_ec   <-  read.delim("https://rest.kegg.jp/link/glycan/ec", header = F, quote = "")


#ID Aliases
cpd_labs      <-  read.delim("https://rest.kegg.jp/list/cpd",     header = F, quote = "")
drug_labs     <-  read.delim("https://rest.kegg.jp/list/drug",    header = F, quote = "")
ec_labs       <-  read.delim("https://rest.kegg.jp/list/ec",      header = F, quote = "")
ko_labs       <-  read.delim("https://rest.kegg.jp/list/ko",      header = F, quote = "")
glycan_labs   <-  read.delim("https://rest.kegg.jp/list/glycan",      header = F, quote = "")


colnames(ko_ec)        = c("ec",  "ko" )
colnames(cpd_drug)     = c("dr",  "cpd")
colnames(cpd_glycan)   = c("gl",  "cpd")

colnames(cpd_ec)       = c("ec",  "cpd")
colnames(drug_ko)      = c("ko",  "dr")
colnames(glycan_ec)    = c("ec", "gl")

colnames(cpd_labs)     = c("cpd", "cpd_labs")
colnames(drug_labs)    = c("dr",  "dr_labs")
colnames(ec_labs)      = c("ec",  "ec_labs")
colnames(ko_labs)      = c("ko",  "ko_labs")
colnames(glycan_labs)  = c("gl",  "gl_labs")

save(ko_ec, cpd_drug, cpd_glycan, 
     cpd_ec, drug_ko, glycan_ec, 
     cpd_labs, drug_labs, ec_labs, ko_labs, glycan_labs, 
     file = "data/ariadne_KEGG.RData", compress = "xz", compression_level = 9)

#load("data/ariadne.RData")



#cross db
cpd_chebi   <-     read.delim("https://rest.kegg.jp/conv/cpd/chebi",      header = F, quote = "")
cpd_pubchem <-     read.delim("https://rest.kegg.jp/conv/cpd/pubchem",    header = F, quote = "")

drug_chebi   <-    read.delim("https://rest.kegg.jp/conv/drug/chebi",      header = F, quote = "")
drug_pubchem <-    read.delim("https://rest.kegg.jp/conv/drug/pubchem",    header = F, quote = "")

glycan_chebi   <-  read.delim("https://rest.kegg.jp/conv/glycan/chebi",      header = F, quote = "")
glycan_pubchem <-  read.delim("https://rest.kegg.jp/conv/glycan/pubchem",    header = F, quote = "")


save(cpd_chebi, cpd_pubchem,  
     drug_chebi, drug_pubchem,
     glycan_chebi, glycan_pubchem, 
     file = "data/ariadne_KEGG_other_DBs.RData", compress = "xz", compression_level = 9)

