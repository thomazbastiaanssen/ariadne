
#library(mia)

#data("Tengeler2020", package = "mia")

#tse <- Tengeler2020

#res <- gsub(".*\\|", "", getFullTaxonomyLabels(rowData(tse)))
#res <- gsub("[a-z]__", "", res)

#graph <- ariadne()

#tax2gmm <- weaveComplexModules(graph, taxname ~ gmm, k = 3, init = res[1:20])

#gmm2tax <- weaveComplexModules(graph, gmm ~ taxname, k = 2)

#exp1 <- weaveComplexModules(graph, taxname ~ gmm, k = 3, init = res[1:5])


#' @export
#' @importFrom Matrix Matrix crossprod colSums
#' @importFrom igraph as_data_frame
weaveComplexModules <- function(graph, by, k, init = NULL, mode = "presence",
    threshold = 1, output.format = "long", verbose = TRUE, ...){
    # Check args
    mode <- match.arg(mode, c("presence", "coverage"))
    output.format <- match.arg(output.format, c("long", "wide"))
    stopifnot(
        "'threshold' must be between 0 and 1." = threshold > 0 && threshold <= 1
    )
    # Extract formula vars
    by.vars <- all.vars(by)
    # Define possible module names
    mod.names <- c("gbm", "gmm")
    # Identify module name
    is.mod <- by.vars %in% mod.names
    # Check that exactly one module name is specified
    if( sum(is.mod) != 1L ){
        stop("Either source or target variable must be a module name.",
            call. = FALSE)
    }
    
    var.idx <- c(mod = which(is.mod), orig = which(!is.mod))
    
    mod.name <- by.vars[var.idx[["mod"]]]
    orig.name <- by.vars[var.idx[["orig"]]]
    
    feat.name <- switch(mod.name, gmm = "ko", gbm = "ko+eggnog+tigr")
    
    edge_df <- as_data_frame(graph, what = "edges")
    url <- unique(edge_df$path[edge_df$from == mod.name])
    
    x <- readLines(url)
    linkmaps <- .process_complex_modules(x)

    if( var.idx[["mod"]] == 1 ){
        init <- unique(linkmaps[["complex2feature"]][["feature"]])
    }
    
    inner.by <- c(feat.name, orig.name)[var.idx] |>
        paste(collapse = "~") |>
        as.formula()
    
    feature2orig <- weavePath(
        graph, inner.by, k, init = init, verbose = verbose, ...
    )
    
    feature2orig <- feature2orig[ , var.idx]
    colnames(feature2orig) <- c("feature", "orig")
    linkmaps[["feature2orig"]] <- feature2orig
    
    mf <- MultiFactor(linkmaps)
    # Print step
    if( verbose ) message(feat.name, " -(GM)-> ", mod.name)
    
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
    
    # Assess module coverage
    out <- Matrix(
        crossprod(m2cp, component2x != 0L) / Matrix::colSums(m2cp),
        sparse = TRUE
    )
    
    if( mode == "presence" ){
        # Compute presence
        out <- out >= threshold
    }
    
    rownames(out) <- levels(mf[["module2component"]][[1L]])
    colnames(out) <- levels(mf[["feature2orig"]][[2L]])
    
    out <- out |>
        as.matrix() |>
        t()
    
    if( output.format == "long" ){
    
        idx <- which(out > 0, arr.ind = TRUE)
        
        out <- data.frame(
            x = rownames(out)[idx[ , 1]],
            y = colnames(out)[idx[ , 2]],
            row.names = NULL
        )
        
        colnames(out) <- by.vars[var.idx]
    }
    return(out)
}


.process_complex_modules <- function(x, br = "///", AND = ",", OR = "\t"){
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
    modules <- list(
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
    return(modules)
}
