mapFiles <- list(
    KO = KO2UniRef90,
    #EF = ...
)

mapTaxonomy <- function(modules, from = "KO"){
  
  map <- mapFiles[[from]]
  
  keep <- names(map) %in% unique(unlist(modules))
  
  map <- map[keep]
  
  # Reduce size for testing
  map <- map[seq(1, 10)]
  map <- lapply(map, function(x) x[seq(1, 10)])
  
  tax.list <- list()
  
  for( i in seq_along(map) ){
    
    module <- map[[i]]
    name <- names(module)[i]
    
    module <- mapUniProt(
      from = "UniProtKB_AC-ID",
      to = "UniRef90",
      query = module,
      columns = c("id", "members")
    )
    
    uniref_ids <- unlist(strsplit(module$Cluster.members, "; "))
    
    tax <- queryUniProt(
      query = uniref_ids,
      fields = c("id", "organism_name", "lineage"),
      collapse = " OR "
    )
    
    tax.list[name] <- list(tax)
    
  }
  
  df <- bind_rows(tax.list, .id = "source")
  
  pattern <- paste(TAXONOMY_RANKS[-2], collapse = "|")
  pattern <- paste0("\\((", pattern, ")\\)")
  
  df <- df[df$Taxonomic.lineage != "", ]
  
  splitted_taxa <- strsplit(df$Taxonomic.lineage, split = ", ")
  
  df$Taxonomic.lineage <- lapply(splitted_taxa, function(x)
    paste0(x[str_detect(x, pattern)], collapse = ", "))
  
  df <- df |>
    separate_wider_delim(cols = "Taxonomic.lineage", delim = ", ",
                         names = TAXONOMY_RANKS[-2], too_few = "align_start")
  
  df[ , TAXONOMY_RANKS[-2]] <- apply(df[ , TAXONOMY_RANKS[-2]], 2L, function(col)
    gsub(" \\(\\w+\\)", "", col))
  
  # This shouldn't be needed
  df[is.na(df)] <- ""
  
  return(df)
}