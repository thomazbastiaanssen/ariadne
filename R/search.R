#' Search paths between resources
#' 
#' @name searchPath
#' 
#' @description
#' \code{searchPath} allows to search the first k-th shortest paths between
#' resources.
#' 
#' @param graph An igraph object.
#' 
#' @param by A formula specifying the path to search.
#' 
#' @param k \code{Numeric scalar}. The kth shortest paths to search.
#'   (Default: \code{1})
#' 
#' @param include \code{Character vector}. Nodes to cross in the path.
#'   (Default: \code{NULL})
#' 
#' @param exclude \code{Character vector}. Nodes to avoid in the path.
#'   (Default: \code{NULL})
#' 
#' @param res.name \code{Character vector}. Names of resources to include in
#'   the graph. (Default: \code{NULL})
#' 
#' @param ... Unused.
#' 
#' @returns
#' A message.
#' 
#' @examples
#' # Retrieve resource graph
#' graph <- ariadne()
#' 
#' # Search first 5 paths from ko to ec
#' searchPath(graph, ko ~ ec, k = 5)
#' 
#' # Search first path including uniref90
#' searchPath(graph, taxname ~ ko, include = "uniref90")
#' 
#' # Search first 5 paths excluding uniref50 and uniref100
#' searchPath(graph, taxname ~ ko, k = 5, exclude = "uniref50")
#' 
#' # Search path for a subset of resources
#' searchPath(graph, uniref90 ~ eggnog, k = 3, res.name = "ChocoPhlAn")
NULL


#' @export
#' @rdname searchPath
#' @importFrom igraph E V k_shortest_paths
setMethod("searchPath", signature = c(graph = "igraph"),
    function(graph, by, k = 1, include = NULL, exclude = NULL, res.name = NULL){
    # Initialise message
    msg <- c()
    # Print paths up to k
    for( j in seq_len(k) ){
        # Add path number
        msg <- c(msg, "Path ", j, ":\n")
        # Get path
        path_df <- .draw_path(graph, by, j, include, exclude, res.name)
        # Add path string
        path_str <- paste0(
            " -(", path_df$source, ")-> ", path_df$to, collapse = ""
        )
        # Add origin
        path_str <- paste0(path_df$from[1], path_str)
        # Add new line 
        msg <- c(msg, path_str, "\n\n")
    }
    # Send message
    message(msg)
    invisible(NULL)
})
