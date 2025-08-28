#' @export
#' @importFrom UniProt.ws mapUniProt queryUniProt
#' @importFrom mia getTaxonomyRanks
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_detect
mapModules <- function(modules, map){
    
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
    
    tax.df <- bind_rows(tax, .id = "source")
    
    tax.ranks <- getTaxonomyRanks()[-2]
    
    pattern <- paste(tax.ranks, collapse = "|")
    pattern <- paste0("\\((", pattern, ")\\)")
    
    tax.df <- tax.df[tax.df$Taxonomic.lineage != "", ]
    
    splitted.tax <- strsplit(tax.df$Taxonomic.lineage, split = ", ")
    
    tax.df$Taxonomic.lineage <- lapply(splitted.tax, function(x)
        paste0(x[str_detect(x, pattern)], collapse = ", "))
    
    tax.df <- tax.df |>
        separate_wider_delim(
            cols = "Taxonomic.lineage", delim = ", ",
            names = tax.ranks, too_few = "align_start"
        )
    
    tax.df[ , tax.ranks] <- apply(tax.df[ , tax.ranks], 2L,
        function(col) gsub(" \\(\\w+\\)", "", col))
    
    # This shouldn't be needed
    tax.df[is.na(tax.df)] <- ""
    
    tax.lab <- mia:::.tax_table2label(tax.df[ , tax.ranks])

    sig.list <- list()
    
    for( i in seq_along(modules) ){
      
        module <- names(modules[i])
        keep <- tax.df$source %in% modules[[i]]
    
        if( any(keep) ){
            sig.list[[module]] <- tax.lab[keep]
        }
    }
    
    return(sig.list)
}
