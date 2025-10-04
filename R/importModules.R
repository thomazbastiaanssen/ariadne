#' Import modules from a file or a database
#'
#' @name importModules
#' @rdname importModules
#'
#' @description
#' \code{importModules} retrieves modules information from a file or database.
#'
#' @param x \code{Character vector}. Path to custom database of module
#'   files or one of the available databases
#'   (\code{c("GBM", "GMM")}).
#' @param br,AND,OR \code{Character scalar}. Used when parsing custom files.
#'  Which characters should be interpreted as the end of a module as well as AND
#'  and OR declarations. (Defaults: \code{'///',',','\t'}.
#'
#' @param verbose \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#'
#' @examples
#'
#' # Import GMM modules
#' modules1 <- importModules("GMM")
#'
#' # Import local module file
#' # modules2 <- importModules("path/to/file")
NULL

ModuleDatabases <- list(
  GMM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GMMs.v1.07.txt",
  GBM = "https://github.com/omixer/omixer-rpmR/raw/refs/heads/main/inst/extdata/GBMs.v1.0.txt"
)


S7::method(importModules, S7::class_character) <- function(
    x, br = "///", AND = ",", OR = "\t",
    module_id = "module", feature_id = "feature", verbose = TRUE
    ){

  if( !is.logical(verbose) ){
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  # Whether to use package or custom modules
  if( x %in% names(ModuleDatabases) ){
    # Cache database
    module_id  <- x
    if(x %in% c("GBM", "GMM")) feature_id <- "ko"
    x <- .getCache(ModuleDatabases[[x]])
  }
  modules <- .import_modules(x, br, AND, OR)

  # Double gsub to set both module and feature names
  names(modules) <- gsub("module", module_id,
                         gsub("feature", feature_id, names(modules)))

  modules <- lapply(
    modules, function(linkmap) `names<-`(
      linkmap, gsub("module", module_id,
              gsub("feature", feature_id, names(linkmap))))
    )
  MultiFactor::MultiFactor(modules)
}

.import_modules <- function(x, br = "///", AND = ",", OR = "\t"){

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

  out <- list(
    module2module_component = data.frame(
      module,
      module_component
      ),
    component2feature_complex = data.frame(
      module_component = rep(module_component, lengths(feature_list)),
      feature_complex
    ),
    feature_complex2feature = data.frame(
      feature_complex = rep(feature_complex, lengths(feature)),
      feature = unlist(feature, use.names = FALSE)
    )
    )

  out

}
