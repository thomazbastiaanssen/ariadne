
mapFiles <- list(
    KO = importMapping("data/KO_uniref90_dict.txt")
    #EF = ...
)

#' @export
#' @importFrom UniProt.ws mapUniProt queryUniProt
#' @importFrom mia getTaxonomyRanks
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_detect
mapTax <- function(modules, from = "KO"){
  
    map <- mapFiles[[from]]
  
    keep <- names(map) %in% unique(unlist(modules))
  
    map <- map[keep]
  
    # Reduce size for testing
    map <- map[seq(1, 10)]
    map <- lapply(map, function(x) x[seq(1, 10)])
  
    tax <- list()
  
    for( i in seq_along(map) ){
    
        item <- map[[i]]
        key <- names(map)[[i]]

        item.members <- mapUniProt(
            from = "UniProtKB_AC-ID",
            to = "UniRef90",
            query = item,
            columns = c("id", "members")
        )

        uniref.ids <- unlist(strsplit(item.members$Cluster.members, "; "))
    
        item.tax <- queryUniProt(
            query = uniref.ids,
            fields = c("id", "organism_name", "lineage"),
            collapse = " OR "
        )
    
        tax[key] <- list(item.tax)
    
    }
  
    df <- bind_rows(tax, .id = "source")
    
    tax.ranks <- getTaxonomyRanks()[-2]
  
    pattern <- paste(tax.ranks, collapse = "|")
    pattern <- paste0("\\((", pattern, ")\\)")
  
    df <- df[df$Taxonomic.lineage != "", ]
  
    splitted_taxa <- strsplit(df$Taxonomic.lineage, split = ", ")
  
    df$Taxonomic.lineage <- lapply(splitted_taxa, function(x)
        paste0(x[str_detect(x, pattern)], collapse = ", "))
  
    df <- df |>
        separate_wider_delim(
            cols = "Taxonomic.lineage", delim = ", ",
            names = tax.ranks, too_few = "align_start"
        )
  
    df[ , tax.ranks] <- apply(df[ , tax.ranks], 2L,
        function(col) gsub(" \\(\\w+\\)", "", col))
  
    # This shouldn't be needed
    df[is.na(df)] <- ""
  
    return(df)
}
