#' Map modules to taxa
#' 
#' \code{mapModules} returns a list of modules containing the taxa that are
#' members of each module. Taxa are derived from uniref ids by querying UniProt
#' SPARQL. Membership is based on whether a taxon meets the criteria specified
#' by a module.
#'
#' @param modules \code{Named character list}. Named list of vectors, where each
#'   vector is a module and its elements are the module members.
#'
#' @param map \code{Named character list}. Named list of vectors, where each
#'   vector is a mapping key and its elements are the mapped values.
#' 
#' @param mode \code{Character scalar}. Specifies the type of modules to expect
#'   as input. It can be one of \code{c("single", "andor")}.
#'   (Default: \code{"single"}).
#' 
#' @param remove.empty \code{Logical scalar}. Should modules with no matching
#'   taxa be removed. (Default: \code{TRUE}).
#' 
#' @param verbose \code{Logical scalar}. Should information on execution be
#'   printed in the console. (Default: \code{TRUE}).
#' 
#' @details
#' The input modules of \code{mapModules} can be either single elements or
#' and/or relationships, depending on whether \code{mode} is \code{"single"} or
#' \code{"andor"}.
#' 
#' @return
#' \code{mapModules} returns a named list of vectors, where each vector is a
#' module and its elements are the members of that module.
#'
#' @examples
#' # Import GBM
#' gbm <- importModules("GBM")
#' 
#' # Import ko-to-uniref90 mapping
#' map <- importMapping("ChocoPhlAn", from = "ko", to = "uniref90")
#' 
#' # Map modules to UniRef90
#' uniref.modules <- mapModules(gbm[seq(3)], map)
#' 
#' # Map modules to taxa
#' tax.modules <- mapModules(gbm[seq(3)], map, mode = "taxonomy")
#' 
#' @name mapModules
NULL

#' @export
#' @rdname mapModules
#' @importFrom BiocParallel bplapply
setMethod("mapModules", signature = c(modules = "list"),
    function(modules, map, mode = "single", remove.empty = TRUE, verbose = TRUE){
        # Check arguments
        if( !is.vector(modules) ){
            stop("'modules' must be a character vector or list of character ",
                "vectors, where each vector corresponds to a module.",
                call. = FALSE)
        }
        if( !is.vector(map) ){
            stop("'map' must be a character vector or list of character ",
                "vectors, where each vector corresponds to a mapping.",
                call. = FALSE)
        }
        if( !mode %in% c("single", "andor") ){
            stop("'mode' must be either single or andor.", call. = FALSE)
        }
        if( !is.logical(remove.empty) ){
            stop("'remove.empty' should be TRUE or FALSE.", call. = FALSE)
        }
        # If a sequence of mappings is provided
        if( all(is.list(unlist(map, use.names = FALSE, recursive = FALSE))) ){
            # Keep only relevant bindings
            keep <- names(map[[1]]) %in% unlist(modules, use.names = FALSE)
            map[[1]] <- map[[1]][keep]
            # Perform mapping recursively
            map <- Reduce(.single_mapping, map)
        }
        # Select mapping method based on module type
        map.method <- switch(mode,
            single = .single_mapping,
            andor = .andor_mapping
        )
        # Map modules and store in modules list
        sig.list <- map.method(modules, map)
        # Remove empty modules
        if( remove.empty ){
            sig.list <- Filter(function(sig) length(sig) > 0, sig.list)
        }
        if( verbose ){
            message(length(unlist(sig.list)), " items were mapped to ",
                length(sig.list), " modules.")
        }
        return(sig.list)
    }
)

# Perform one-to-one mapping
.single_mapping <- function(x, y){
    # Filter y keys that match x values
    y <- y[names(y) %in% unique(unlist(x, use.names = FALSE))]
    if( length(y) == 0 ){
        stop("Items in 'x' did not match any item in 'y'.", call. = FALSE)
    }
    # Store taxa in modules list
    z <- bplapply(x, function(values)
        unlist(y[names(y) %in% values], use.names = FALSE))
    return(z)
}

# Perform and/or mapping (reaction pathway modules)
.andor_mapping <- function(x, y){
    # Filter y keys that match x values
    y <- y[names(y) %in% unique(unlist(x, use.names = FALSE))]
    # Find functions for each taxon
    linkmap <- as.linkmap(y)
    tax <- split(linkmap$x, linkmap$y)
    # Store taxa in modules list
    z <- lapply(x, function(module) {
        members <- vapply(tax, function(tax.item) {
            all(vapply(module, function(comp) any(comp %in% tax.item), logical(1L)))
        }, logical(1L))
        names(tax)[members]
    })
    return(z)
}

# .andor2regex <- function(modules){
#     modules <- lapply(gbm, function(module){
#         module <- vapply(module, function(item){
#             paste(vapply(item, function(comp){
#                 if (length(comp) > 1L) comp <- sprintf("(?=.*?%s)", comp)
#                 comp <- paste(comp, collapse = "")
#                 return(comp)
#             }, character(1L)),
#             collapse = "|")
#         },
#         character(1L))
#         module <- paste(sprintf("(?=.*?%s)", module), collapse = "")
#         return(module)
#     })
#     return(modules)
# }
