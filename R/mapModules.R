#' Map modules to taxa
#' @name mapModules
#' @rdname mapModules
#'
#' @description
#' \code{mapModules} returns a list of modules containing the taxa that are
#' members of each module. Taxa are derived from uniref ids by querying UniProt
#' SPARQL. Membership is based on whether a taxon meets the criteria specified
#' by a module.
#'
#' @param map \code{MultiFactor}. Typically produced by `importModules()`.
#'
#' @param x \code{matrix}, or object that can be coerced to `matrix`, such as a
#'     `data.frame`, with features as rows and samples as columns.
#'
#' @param coverage.threshold \code{Numeric scalar}. Minimum proportion of
#'     components (i.e., reaction steps) required to be considered present.
#'     bounded between 0-1. (Default: 0.8)
#' @param method \code{Character vector}.
#' @param verbose \code{Logical scalar}. Should information on execution be
#'     printed in the console. (Default: \code{TRUE}).
#'
#' @returns \code{mapModules} returns a numeric matrix with modules as rows and
#'     samples as columns, with content depending on desired `mode` argument.
#'
#' @examples
#' # Import GBM
#' map <- importModules("GBM")
#'
#' # sparse feature table with 20 random samples
#' x_present <- replicate(20,
#'     sample(
#'     c(TRUE, FALSE), size = nlevels(map$ko_complex2ko[[2L]]), replace = TRUE
#'     ),
#'     simplify = TRUE
#' )
#' x <- x_present * matrix(
#'     rnorm(n = prod(dim(x_present))),
#'     nrow = NROW(x_present), ncol = 20
#' )^2
#'
#' # Map modules to feature table
#' mapModules(map, x, method = "coverage")
#' mapModules(map, x, method = "sum")
#' mapModules(map, x, method = "count")
#' mapModules(map, x, method = "presence")
#'
NULL


#' @importFrom BiocParallel bplapply
#' @importFrom Matrix crossprod
#' @importFrom MultiFactor weave
#' @noRd
S7::method(mapModules, MultiFactor) <- function(
        map, x, method = c("sum", "count", "coverage", "presence"),
        verbose = TRUE, coverage.threshold = 0.8
        ) {

    # Check arguments
    if( !is.logical(verbose) ){
        stop("'verbose' should be TRUE or FALSE.", call. = FALSE)
    }
    if(coverage.threshold <= 0L | coverage.threshold > 1L ) {
        stop("'coverage.threshold' should be between 0-1.", call. = FALSE)
    }
    if(!inherits(x, "matrix")) {
        x <- as.matrix(x)
    }

    method <- match.arg(arg = method, c("sum", "count", "coverage", "presence"))

    c2f_ind <- grep("_complex2", names(map))
    c2c_ind <- grep("^component2.*_complex$", names(map))
    m2c_ind <- grep("_component$", names(map))

    # Check for presence of all complex components
    complex2x   <- .map_AND_complex(map[[c2f_ind]], x)
    c2c         <- as.matrix(map[[c2c_ind]], terms = c(2L, 1L))
    component2x <- Matrix::crossprod(c2c, complex2x, boolArith = TRUE)
    # Assess module coverage
    module_coverage <- .map_coverage(map[[m2c_ind]], component2x)

    if( verbose ){
        message(NROW(x), " features were mapped to ",
                NROW(module_coverage), " modules.")
    }
    if( method == "coverage" ) {
        rownames(module_coverage) <- levels(map[[m2c_ind]][[1L]])
        return(module_coverage)
    }
    module_present <- module_coverage >= coverage.threshold
    if( method == "presence" ) {
        rownames(module_present) <- levels(map[[m2c_ind]][[1L]])
        return(module_present)
    }

    module_name <- names(map[[m2c_ind]])[[1L]]
    x_name <- names(map[[c2f_ind]])[[2L]]

    mod2x <- as.matrix(MultiFactor::weave(map, .by = c(x_name, module_name)))

    if( method == "sum" ) {
        sum_table <- Matrix::crossprod(mod2x, x) * module_present

        rownames(sum_table) <- levels(map[[m2c_ind]][[1L]])
        return(sum_table)
    }

    if( method == "count" ) {
        count_table <- Matrix::crossprod(mod2x, x != 0L) * module_present

        rownames(count_table) <- levels(map[[m2c_ind]][[1L]])
        return(count_table)
    }

}

# x = module limkmap, with complex on the left (1L), features on the right (2L).
# y = feature table
# return = bool matrix indicating full complex presence, row, per sample, col.
#
.map_AND_complex <- function(x, y) {
    x_mat <- as.matrix(x, terms = c(2L, 1L))
    Matrix::Matrix(
        Matrix::crossprod(x_mat, y!=0L) >= Matrix::colSums(x_mat),
        sparse = TRUE
    )
}

# x = module linkmap, with module on the left (1L), component on the right (2L).
# y = feature table
# return = numeric matrix indicating module coverage as prop, row, per sample, col.
#
.map_coverage <- function(x, y) {
    x_mat <- as.matrix(x, terms = c(2L, 1L))
    Matrix::Matrix(
        Matrix::crossprod(x_mat, y!=0L) / Matrix::colSums(x_mat),
        sparse = TRUE
    )
}


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

#' @importFrom MultiFactor as.LinkMap
# Perform and/or mapping (reaction pathway modules)
.andor_mapping <- function(x, y){
    # Filter y keys that match x values
    y <- y[names(y) %in% unique(unlist(x, use.names = FALSE))]
    # Find functions for each taxon
    linkmap <- as.LinkMap(y)
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



# Query taxonomies based on uniref ids from UniProt using SPARQL
# x = Character vector of uniref IDs
.querySPARQL <- function(x, graph = rdflib$Graph()){
    # Collapse UniRef90 ids into long string
    uniref.ids <- paste0("uniref:", x, collapse = " ")
    # Define first part of query
    query_part1 <- "
        PREFIX uniprot: <http://purl.uniprot.org/core/>
        PREFIX uniref: <http://purl.uniprot.org/uniref/>
        # Outer query to get final names
        SELECT ?name
        WHERE {
            SERVICE <https://sparql.uniprot.org/> {
            {   # Inner query to get distinct taxa
                SELECT DISTINCT ?taxId
                WHERE {
                    # List UniRef90 ids
                    VALUES ?unirefId {
        "
    # Define second part of query
    query_part2 <- "
                    }
                    # Bind UniRef90 ids to cluster members
                    ?unirefId uniprot:member ?member.
                    # Bind cluster members to taxa
                    ?member uniprot:organism ?taxId.
                }
            }
            # Bind taxa to scientific names and ranks
            ?taxId uniprot:scientificName ?sciName; uniprot:rank ?rank.
            # Process name to prefix__taxon format
            BIND(lcase(substr(strafter(str(?rank), 'Rank_'), 1, 1)) as ?prefix)
            BIND(concat(?prefix, '__', ?sciName) as ?name)
            }
        }
        "
    # Build query
    query <- paste0(query_part1, "\t\t\t", uniref.ids, query_part2)
    # Execute query
    qres <- graph$query(query)
    # Convert name bindings to vector
    tax.vec <- vapply(qres$bindings, function(binding)
        binding[["name"]], character(1L))
    return(tax.vec)
}
