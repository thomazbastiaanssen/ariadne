#' Add resources to graph
#'
#' @name addResource
#'
#' @description
#' \code{addResource} includes a user-defined resource to the ariadne graph.
#' 
#' @param graph An igraph object.
#' 
#' @param file \code{Character scalar}. The path to a local or remote file that
#'   stores the linkmap between two features. It must be a table with two named
#'   columns (features) and each row representing a binding between features.
#' 
#' @param res.name \code{Character scalar}. The name of the resource that will
#'   appear in the output graph. (Default: \code{"Custom"})
#' 
#' @param force \code{Logical scalar}. Whether \code{file} should overwrite any
#'   previously cached files with the same name. (Default: \code{FALSE})
#' 
#' @param ... Additional arguments passed to \code{\link[data.table:fread]{fread}}.
#' 
#' @returns An igraph object.
#' 
#' @examples
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Set URL to custom resource
#' url <- "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-12/"
#' url <- paste0(url, "genotype.csv.gz")
#' 
#' # Add resource to ariadne graph
#' graph <- addResource(
#'     graph,
#'     url,
#'     res.name = "AMR",
#'     select = c("taxon_id", "antibiotic_ontology"),
#'     col.names = c("taxid", "aro")
#' )
#' 
#' # Search for newly added path
#' searchPath(graph, taxid ~ aro)
#' 
#' # Plot first 5 paths excluding uniref50 and uniref100
#' plotPath(graph, taxid ~ aro, focus = TRUE)
#' 
#' # Find antibiotics related to E. coli (NCBI 562)
#' tax2aro <- weavePath(graph, taxid ~ aro, init = 562)
NULL

#' @export
#' @rdname addResource
#' @importFrom igraph as_data_frame graph_from_data_frame
#' @importFrom BiocFileCache bfcquery
#' @importFrom dplyr bind_rows
#' @importFrom data.table fread
setMethod("addResource", signature = c(graph = "igraph"),
    function(graph, file, res.name = "Custom", force = FALSE, ...){
    # Break graph into edges and nodes
    graph_df <- as_data_frame(graph, what = "both")
    edge_df <- graph_df$edges
    node_df <- graph_df$vertices
    # Initialise cache
    bfc <- .init_cache()
    # Build resource name
    rname <- file.path(res.name, basename(file))
    # Check if preprocessed file is cached
    cached <- bfcquery(bfc, rname)$rpath
    # Return preprocessed file if available
    if( !force && length(cached) > 0L ){
        stop("This link for ", res.name, " already exists. Set 'force' to TRUE",
            " to overwrite it.", call. = FALSE)
    }
    # Import linkmap
    linkmap <- fread(file, ...)
    # Check that linkmap has two columns
    if( ncol(linkmap) != 2L ){
        stop("'file' must point to a two-column linkmap.", call. = FALSE)
    }
    # Retrieve feature names from colnames
    cus_vars <- colnames(linkmap)
    # Find if some variables are new
    new_vars <- setdiff(cus_vars, node_df$name)
    # Check that at least one feature is known
    if( length(new_vars) == 2L ){
        stop("At least one feature must be in 'graph'.", call. = FALSE)
    }
    # Add new edge for user-defined linkmap
    edge_df <- rbind(edge_df, data.frame(
        from = cus_vars[1], to = cus_vars[2], source = res.name, url = file
    ))
    # Find duplicate edges (originally none)
    is_dup <- edge_df |>
        .get_edge_keys() |>
        duplicated(fromLast = TRUE)
    # Check that new edge is not duplicate
    if( !force && any(is_dup) ){
        stop("This link for ", res.name, " already exists. Set 'force' to TRUE",
            " to overwrite it.", call. = FALSE)
    }
    # Update duplicate with user-defined binding
    edge_df <- edge_df[!is_dup, ]
    # Add nodes for new variables (if any)
    node_df <- bind_rows(node_df, data.frame(name = new_vars))
    # If resource is new
    if( is.null(node_df[[res.name]]) ){
        # Create column for resource-specific names
        node_df[[res.name]] <- NA
    }
    # Add resource-specific names
    node_df[node_df$name %in% cus_vars, res.name] <- cus_vars
    # Store linkmap in cache as parquet file
    .add2cache(linkmap, rname, bfc)
    # Combine edges and nodes data into graph
    graph <- graph_from_data_frame(edge_df, vertices = node_df)
    return(graph)
})
