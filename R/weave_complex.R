
#' @export
#' @rdname weavePath
#' @importFrom stats as.formula
setMethod("weaveComplex", signature = c(graph = "data.frame"),
    function(graph, init = NULL, prune = TRUE, use.names = TRUE,
    threshold = NULL, verbose = TRUE, timeout = 1e6, ...){
    
    by <- c(graph$from[1L], graph$to[nrow(graph)]) |>
        paste(collapse = "~") |>
        as.formula()
    
    graph <- .graph_from_path_df(graph)
    
    linkmap <- weaveComplex(
        graph, by, init = init, prune = prune, use.names = use.names,
        threshold = threshold, verbose = verbose, timeout = timeout, ...
    )
    return(linkmap)
})


#' @export
#' @rdname weavePath
#' @importFrom igraph as_data_frame
#' @importFrom stats as.formula reformulate
#' @importFrom Matrix summary
setMethod("weaveComplex", signature = c(graph = "igraph"),
    function(graph, by, k = 1, include = NULL, exclude = NULL, res.name = NULL,
    init = NULL, prune = TRUE, use.names = TRUE, threshold = NULL,
    verbose = TRUE, timeout = 1e6, ...){
    # Check threshold
    if( !is.null(threshold) && (!is.numeric(threshold) ||
        length(threshold) != 1L || threshold <= 0 || threshold > 1) ){
        stop("'threshold' must be a number between 0 and 1.", call. = FALSE)
    }
    # Extract formula vars
    by.vars <- all.vars(by)
    # Identify module name
    origin <- by.vars[1L]
    target <- by.vars[2L]
    # Define complex modules
    complex_modules <- c("gbm", "gmm")
    
    if( target %in% complex_modules ){
    
        edge_df <- as_data_frame(graph, what = "edges")
        
        inter_name <- edge_df$to[edge_df$from == target]
        url <- edge_df$url[edge_df$from == target]
        
        linkmaps <- .process_complex_modules(url, output.format = "list")
        
        inner_by <- reformulate(inter_name, origin)
        
    }else{
        inner_by <- by
    }
    # Build MultiFactor from path linkmaps
    mf <- .build_path_mf(
        graph, inner_by, k, include, exclude, res.name,
        init, prune, TRUE, verbose, timeout, ...
    )
    
    if( target %in% complex_modules ){
        # Print step
        if( verbose ) message(inter_name, " -(GM)-> ", target)
        # Weave first part of the path
        origin2feature <- weave(mf, inner_by) |> as.data.frame()
        # Add names to linkmap columns
        colnames(origin2feature) <- c("origin", "feature")
        # Include origin2feature linkmap in list
        linkmaps[["origin2feature"]] <- origin2feature
        # Construct MultiFactor from linkmaps
        mf <- MultiFactor(linkmaps)
        # Map features to complex modules
        mat <- .map_complex_modules(mf)
    }else{
        
        outer_by <- colnames(mf)[c(1L, ncol(mf) - 1)] |>
            paste(collapse = "~") |>
            as.formula()
        
        feature2origin <- mf |>
            weave(outer_by) |>
            as.matrix(terms = c(2L, 1L))
        
        feature2module <- as.matrix(mf[[nrow(mf)]])
        
        module2origin <- crossprod(feature2module, feature2origin != 0L)
        
        mat <- module2origin / colSums(feature2module)
    }
    # Convert to matrix object
    out <- summary(mat)
    # If defined, subset values above threshold
    if( !is.null(threshold) ) out <- out[out$x >= threshold, ]
    # Convert to linkmap
    out <- data.frame(
        x = as.factor(colnames(mat)[out$j]),
        y = as.factor(rownames(mat)[out$i]),
        z = out$x, row.names = NULL
    )
    # Add colnames
    colnames(out) <- c(origin, target, "cov")
    # Add feature names
    if( use.names ){
        name_links <- linkNames(graph, target, out[[2L]], verbose = verbose)
        out[paste0(target, ".name")] <- as.factor(name_links[[2L]])
    }
    return(out)
})


#' @importFrom Matrix crossprod colSums
.map_complex_modules <- function(mf){
    # Retrieve matrices from MultiFactor
    mats <- lapply(mf, as.matrix, terms = c(2L, 1L))
    # Link binarised origin ids to complexes
    complex2origin <- crossprod(
        mats[["complex2feature"]], mats[["origin2feature"]] != 0L
    )
    # Ensure that a given origin id links to all members of each complex
    complex2origin <- complex2origin == colSums(mats[["complex2feature"]])
    # Link origin ids to components
    component2origin <- crossprod(
        mats[["component2complex"]], complex2origin, boolArith = TRUE
    )
    # Link binarised origin ids to modules
    module2origin <- crossprod(
        mats[["module2component"]], component2origin != 0L
    )
    # Estimate coverage dividing by total of components in each module
    module2origin <- module2origin / colSums(mats[["module2component"]])
    return(module2origin)
}


#' @importFrom readr read_lines
#' @importFrom MultiFactor MultiFactor weave
.process_complex_modules <- function(
    x, br = "///", AND = ",", OR = "\t", output.format = "linkmap", ...){
    # Read file content
    x <- read_lines(x, ...)
    # Identify break lines
    v_br <- x == br
    # Split content by breaks, excluding break lines themselves
    line_content <- split(x[!v_br], cumsum(v_br)[!v_br])
    # Extract keys (first line of each block)
    keys <- vapply(line_content, `[`, 1L, FUN.VALUE = "", USE.NAMES = FALSE)
    # Extract values (all lines except first in each block)
    values <- lapply(line_content, `[`, -1L)
    # Replace tabs with spaces in keys, repeated for each value line
    module <- gsub("\t.*", "", rep(keys, lengths(values, use.names = FALSE)))
    # Create unique module_component identifiers
    module_component <- paste0(
        module, "_part_",
        unlist(lapply(rle(module)$lengths, seq), use.names = FALSE)
    )
    # Split feature list by tab character
    feature_list <- strsplit(
        unlist(values, recursive = TRUE, use.names = FALSE), "\t", fixed = TRUE
    )
    names(feature_list) <- module_component
    # Flatten feature complex list and split by comma to get individual features
    feature_complex <- unlist(feature_list, use.names = FALSE)
    feature <- strsplit(feature_complex, ",", fixed = TRUE)
    # Create and return structured list of data frames
    out <- list(
        module2component = data.frame(module, component = module_component),
        component2complex = data.frame(
            component = rep(module_component, lengths(feature_list)),
            complex = feature_complex
        ),
        complex2feature = data.frame(
            complex = rep(feature_complex, lengths(feature)),
            feature = unlist(feature, use.names = FALSE)
        )
    )
    if( output.format == "linkmap" ){
        out <- out |>
            MultiFactor() |>
            weave(module ~ feature) |>
            data.frame()
    }
    return(out)
}
