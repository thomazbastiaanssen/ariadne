#' Visualise modules network
#' 
#' @name plotModules
#' 
#' @examples
#' graph <- ariadne()
#' ec2gmm <- weaveComplex(graph, ec ~ gmm)
NULL

#' @export
#' @rdname plotModules
#' @importFrom data.table as.data.table melt
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggraph ggraph scale_edge_colour_gradient2 geom_node_text
#' @importFrom ggplot2 scale_colour_manual guides guide_legend
setMethod("plotModules", signature = c(modules = "data.frame"),
    function(modules, edge.type = "diagonal", show.labels = TRUE, ...){
    
    if( is.null(attr(modules, "path.meta")) ){}
    
    modules <- as.data.table(modules)
    if( !"cov" %in% colnames(modules) ) modules$cov <- 1
    
    node_df <- melt(
        modules[ , c(1, 2)], measure.vars = c(1, 2),
        variable.name = "type", value.name = "name",
        variable.factor = TRUE, value.factor = TRUE
    )
    
    node_df <- unique(node_df[ , c(2, 1)])
    graph <- graph_from_data_frame(modules, vertices = node_df)
    
    geom_edge <- eval(parse(text = paste0("ggraph::geom_edge_", edge.type)))
    
    p <- ggraph(graph, ...) +
        geom_edge(aes(colour = .data$cov)) +
        geom_node_point(aes(colour = .data$type)) +
        scale_colour_manual(values = c("darkorange", "#06B4B4")) +
        scale_edge_colour_gradient2(
            low = "white", mid = "grey80", high = "red",
            midpoint = 0.5, limits = c(0, 1)
        ) +
        theme_void() +
        guides(
            colour = guide_legend(order = 1),
            edge_colour = guide_legend(order = 2)
        )
    
    if( show.labels ){
        p <- p + geom_node_text(aes(label = .data$name), vjust = 1.8, size = 2)
    }
    
    return(p)
})

#mods <- go2kegg_path

#edge_df <- attr(mods, "path.meta") |>
#    rbindlist(use.names = FALSE)

#node_df <- mods |>
#    attr("path.meta") |>
#    levels() |>
#    lapply(as.data.table) |>
#    rbindlist(idcol = TRUE)

#node_df <- node_df[ , c(2, 1)]

#colnames(node_df) <- c("name", "type")

#node_df <- node_df[!duplicated(node_df$name), ]

#gr <- graph_from_data_frame(edge_df, vertices = node_df)

#p <- ggraph(gr) +
#    geom_edge_diagonal(colour = "grey80", alpha = 1/10) +
#    geom_node_point(aes(colour = .data$type)) +
#    theme_void() +
#    guides(
#      colour = guide_legend(order = 1),
#      edge_colour = guide_legend(order = 2)
#    )
