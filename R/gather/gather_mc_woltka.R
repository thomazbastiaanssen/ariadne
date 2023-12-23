#Aliases
enzrxn_reaction <- read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/enzrxn-to-reaction.map", header = F, quote = "")
protein_enzrxn  <- read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/protein-to-enzrxn.map",  header = F, quote = "")
protein_gene    <- read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/protein-to-gene.map",    header = F, quote = "")
protein_go      <- read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/protein-to-go.map",      header = F, quote = "")
reaction_ec     <- read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/reaction-to-ec.map",     header = F, quote = "")

#Bridges
reaction_left_compound   <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/reaction-to-left_compound.map",  header = F, quote = "")
reaction_right_compound  <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/reaction-to-right_compound.map", header = F, quote = "")

#ID Aliases
protein_name      <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/protein_name.txt",     header = F, quote = "")
reaction_name     <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/reaction_name.txt",    header = F, quote = "")
compound_inchi    <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/compound_inchi.txt",   header = F, quote = "")
compound_name     <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/compound_name.txt",    header = F, quote = "")
compound_smiles   <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/compound_smiles.txt",  header = F, quote = "")
enzrxn_name       <-  read.delim("http://ftp.microbio.me/pub/wol2/function/metacyc/enzrxn_name.txt",      header = F, quote = "")


save(enzrxn_reaction, protein_enzrxn, protein_gene, protein_go, reaction_ec,
     reaction_left_compound, reaction_right_compound, 
     protein_name, reaction_name, compound_inchi, compound_name, compound_smiles, enzrxn_name, 
     file = "data/ariadne_MC.RData", compress = "xz", compression_level = 9)

#load("data/ariadne_MC.RData")
