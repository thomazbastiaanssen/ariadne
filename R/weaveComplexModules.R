
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
#' @importFrom Matrix Matrix
#' @importFrom igraph as_data_frame
weaveComplexModules <- function(graph, by, k, init = NULL, mode = "presence",
    threshold = 1, output.format = "long", ...){
    # Extract formula vars
    by.vars <- all.vars(by)
    # Define possible module names
    mod.names <- c("gbm", "gmm")
    # Identify module and origin names
    mod.name <- intersect(mod.names, by.vars)
    orig.name <- setdiff(by.vars, mod.names)
    # Check that exactly one module name is specified
    if( length(mod.name) != 1L ){
        stop("Either source or target variable must be a module name.",
            call. = FALSE)
    }
    
    edge_df <- as_data_frame(graph, what = "edges")
    url <- unique(edge_df$path[edge_df$from == mod.name])
    
    x <- readLines(url)
    linkmaps <- .process_complex_modules(x)
    
    feat.name <- switch(mod.name, gmm = "ko", gbm = "ko+eggnog+tigr")

    if( by.vars[1] == mod.name ){
        init <- unique(linkmaps[["complex2feature"]][["feature"]])
        from <- feat.name
        to <- orig.name
    }else{
        from <- orig.name
        to <- feat.name
    }
    
    orig2feature <- weavePath(
        graph, as.formula(paste(from, "~", to)), k, init = init, ...
    )
    
    if( from == feat.name ){
        colnames(orig2feature) <- c("feature", "orig")
    }else{
        colnames(orig2feature) <- c("orig", "feature")
    }
    
    linkmaps[["orig2feature"]] <- orig2feature
    mf <- MultiFactor(linkmaps)
    
    col.order <- c(2L, 1L)
    orig2f <- as.matrix(mf[["orig2feature"]], terms = col.order)
    ct2cx <- as.matrix(mf[["component2complex"]], terms = col.order)
    cx2f <- as.matrix(mf[["complex2feature"]], terms = col.order)
    m2cp <- as.matrix(mf[["module2component"]], terms = col.order)
    
    complex2x <- Matrix(crossprod(cx2f, orig2f != 0L) >= colSums(cx2f), sparse = TRUE)
    component2x <- crossprod(ct2cx, complex2x, boolArith = TRUE)
    
    # Assess module coverage
    out <- Matrix(
        crossprod(m2cp, component2x != 0L) / colSums(m2cp), sparse = TRUE
    )
    
    if( mode == "presence" ){
        # Compute presence
        out <- out >= threshold
    }
    
    rownames(out) <- levels(mf[["module2component"]][[1L]])
    colnames(out) <- levels(mf[["orig2feature"]])[[1L]]
    
    out <- out |>
        t() |>
        as.matrix()
    
    if( output.format == "long" ){
    
        idx <- which(out > 0, arr.ind = TRUE)
        
        out <- data.frame(
            x = rownames(out)[idx[ , 1]],
            y = colnames(out)[idx[ , 2]],
            row.names = NULL
        )
        # Fix from and to in general
        # colnames()
    }
    return(out)
}
