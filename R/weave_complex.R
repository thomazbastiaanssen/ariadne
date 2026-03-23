
#' @importFrom igraph as_data_frame
S7::method(weaveComplex, igraph) <- function(graph, by, k = 1, include = NULL,
    exclude = NULL, init = NULL, prune = TRUE, use.names = TRUE,
    mode = "presence", threshold = 1, verbose = TRUE, timeout = 1e6, ...){
    # Check mode
    mode <- match.arg(mode, c("presence", "coverage"))
    # Check threshold
    if( !is.numeric(threshold) || length(threshold) != 1L ||
        threshold <= 0 || threshold > 1 ){
        stop("'threshold' must be a number between 0 and 1.", call. = FALSE)
    }
    # Extract formula vars
    by.vars <- all.vars(by)
    # Define possible module names
    mod.names <- c("gbm", "gmm")
    # Identify module name
    is.mod <- by.vars %in% mod.names
    # Check that exactly one module name is specified
    if( sum(is.mod) != 1L ){
        stop("Exactly one side of 'by' must be a module name.", call. = FALSE)
    }
    
    var.idx <- c(mod = which(is.mod), orig = which(!is.mod))
    
    mod.name <- by.vars[var.idx[["mod"]]]
    orig.name <- by.vars[var.idx[["orig"]]]
    
    graph_df <- as_data_frame(graph, what = "both")
    edge_df <- graph_df$edges

    feat.name <- edge_df$to[edge_df$from == mod.name]
    url <- edge_df$url[edge_df$from == mod.name]
    
    x <- readLines(url)
    linkmaps <- .process_complex_modules(x, output.format = "list")
    
    if( var.idx[["mod"]] == 1L ){
        init <- unique(linkmaps[["complex2feature"]][["feature"]])
    }
    
    inner_by <- c(feat.name, orig.name)[var.idx] |>
        paste(collapse = "~") |>
        as.formula()
    
    feature2orig <- weavePath(
        graph, inner_by, k, include, exclude,
        init, prune, FALSE, verbose, timeout, ...
    )
    
    feature2orig <- feature2orig[ , var.idx]
    colnames(feature2orig) <- c("feature", "orig")
    linkmaps[["feature2orig"]] <- feature2orig
    # Construct MultiFactor from linkmaps
    mf <- MultiFactor(linkmaps)
    # Print step
    if( verbose ) message(feat.name, " -(GM)-> ", mod.name)
    # Map features to complex modules
    out <- .map_modules(mf)
    # If presence is set
    if( mode == "presence" ){
        # Convert to adjacency matrix
        out <- out >= threshold
    }
    # Set dimnames
    rownames(out) <- levels(mf[["module2component"]][[1L]])
    colnames(out) <- levels(mf[["feature2orig"]][[2L]])
    # Convert to matrix object
    out <- out |>
        as.matrix() |>
        t()
    # Find indices of non-null values
    idx <- which(out > 0, arr.ind = TRUE)
    # Convert to linkmap
    out <- data.frame(
        x = rownames(out)[idx[ , 1L]],
        y = colnames(out)[idx[ , 2L]],
        row.names = NULL
    )
    # Add colnames
    colnames(out) <- c(orig.name, mod.name)
    # Add feature names
    out <- if( use.names ) .id2name(graph_df, out) else out
    return(out)
}


#' @importFrom Matrix Matrix crossprod colSums
.map_modules <- function(mf){
    # Retrieve matrices from MultiFactor
    col.order <- c(2L, 1L)
    orig2f <- as.matrix(mf[["feature2orig"]])
    ct2cx <- as.matrix(mf[["component2complex"]], terms = col.order)
    cx2f <- as.matrix(mf[["complex2feature"]], terms = col.order)
    m2cp <- as.matrix(mf[["module2component"]], terms = col.order)
    
    complex2x <- Matrix(
        crossprod(cx2f, orig2f != 0L) >= Matrix::colSums(cx2f),
        sparse = TRUE
    )
    
    component2x <- crossprod(ct2cx, complex2x, boolArith = TRUE)
    # Compute module coverage
    out <- Matrix(
        crossprod(m2cp, component2x != 0L) / Matrix::colSums(m2cp),
        sparse = TRUE
    )
}


#' @importFrom MultiFactor MultiFactor weave
.process_complex_modules <- function(
    x, br = "///", AND = ",", OR = "\t", output.format = "linkmap"){
    # Identify break lines
    v_br <- x == br
    # Split content by breaks, excluding break lines themselves
    line.content <- split(x[!v_br], cumsum(v_br)[!v_br])
    # Extract keys (first line of each block)
    keys <- vapply(line.content, `[`, 1L, FUN.VALUE = "", USE.NAMES = FALSE)
    # Extract values (all lines except first in each block)
    values <- lapply(line.content, `[`, -1L)
    # Replace tabs with spaces in keys, repeated for each value line
    module <- gsub("\t.*", "", rep(keys, lengths(values, use.names = FALSE)))
    # Create unique module_component identifiers
    module_component <- paste0(
        module, "_part_",
        unlist(lapply(rle(module)$lengths, seq), use.names = FALSE)
    )
    # Split feature list by tab character
    feature_list <- strsplit(
        unlist(values, recursive = TRUE, use.names = FALSE), "\t"
    )
    names(feature_list) <- module_component
    # Flatten feature complex list and split by comma to get individual features
    feature_complex <- unlist(feature_list, use.names = FALSE)
    feature <- strsplit(feature_complex, ",")
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
