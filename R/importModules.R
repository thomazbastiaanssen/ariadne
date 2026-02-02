

.import_modules <- function(x, br = "///", AND = ",", OR = "\t",
                            module_id = "module", feature_id = "feature"){

  # Read the file content
  line.content <- readLines(x)

  v_br <- line.content == br
  line.content <- split(line.content[!v_br], cumsum(v_br)[!v_br])

  keys <- vapply(line.content, `[`, 1L, FUN.VALUE = "", USE.NAMES = FALSE)
  # Extract values
  values <- lapply(line.content, `[`, -1L)

  module <-
    gsub("\t", " ", rep(keys, lengths(values, use.names = FALSE)), fixed = TRUE)

  module_component <- paste0(
    module, "_part_",
    unlist(lapply(rle(module)$lengths, seq), use.names = FALSE)
  )

  feature_list <- strsplit(unlist(values, recursive = TRUE, use.names = FALSE), "\t")
  names(feature_list) <- module_component

  feature_complex <- unlist(feature_list, use.names = FALSE)
  feature <-  strsplit(feature_complex, ",")

  modules <- list(
      module2module_component = data.frame(module, module_component),
      component2feature_complex = data.frame(
          module_component = rep(module_component, lengths(feature_list)),
          feature_complex
      ),
      feature_complex2feature = data.frame(
          feature_complex = rep(feature_complex, lengths(feature)),
          feature = unlist(feature, use.names = FALSE)
      )
  )

  # Double gsub to set both module and feature names
  names(modules) <- gsub(
      "module", module_id, gsub("feature", feature_id, names(modules))
  )

  modules <- lapply(
      modules, function(linkmap) `names<-`(
          linkmap, gsub("module", module_id,
              gsub("feature", feature_id, names(linkmap))))
  )
  
  return(modules)
}
