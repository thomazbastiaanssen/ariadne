
#' @importFrom BiocParallel bpworkers
.get_batches <- function(x, batch.size, workers, factor){
    
    if( is.null(x) ){
        return(list(c(1L, 1L)))
    }
    
    xlen <- length(x)
    
    if( is.null(workers) ){
        workers <- bpworkers()
    }
    
    batch.num <- min(ceiling(xlen / batch.size), factor * workers)
    batch.size <- ceiling(xlen / batch.num)
    
    ranges <- lapply(seq_len(batch.num), function(i) {
        start <- (i - 1) * batch.size + 1
        end <- min(i * batch.size, xlen)
        c(start, end)
    })
    
    return(ranges)
}
