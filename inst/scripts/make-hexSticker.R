
library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(systemfonts)
library(tidygraph)

# Set random seed based on date
seed <- as.numeric(format(Sys.time(), "%y%m%d"))
set.seed(seed)

asset_dir <- "inst/assets/"

# Custom functions ----

hex_graph <- function(size = 10) {
    size <- floor(size)
    # height, width, radius, side (apothem)
    h <- 2 * size
    w <- round(h * sqrt(3) / 2)
    r <- h / 2
    s <- w / 2
    
    x_range <- seq_len(w) - 1L
    y_range <- seq_len(h) - 1L
    
    y <- seq(from = r / 2, to = 0, length.out = size)
    yy <- c(y, rev(y))
    xx <- seq(from = 0,  to = w, length.out = 2 * size)
    x <- xx[seq_len(size)]
    
    corners_drop <- .carve_corners(size, y, x, h, w)
    
    textbox <- expand.grid(
        y = round(c(0, 1) + h / 3),
        x = seq(from = floor(r / 2) + 1, to = ceiling(w - 1 - r / 2) - 1)
    )
    tb_entry <- textbox[range(1, NROW(textbox)), ]
    textbox <- textbox[-range(1, NROW(textbox)), ]
    
    drop_df <- unique(rbind(corners_drop, textbox))
    
    vert_x <- rep(x_range, h)
    vert_y <- rep(y_range, each = w)
    
    drop <- paste(vert_x, vert_y) %in% paste(drop_df$x, drop_df$y)
    left_entry <- vert_x == tb_entry[["x"]][1L] & vert_y == tb_entry[["y"]][1L]
    right_entry <- vert_x == tb_entry[["x"]][2L] & vert_y == tb_entry[["y"]][2L]
    
    g <- make_lattice(c(w, h))
    
    vertex_attr(g, "x") <- vert_x
    vertex_attr(g, "y") <- vert_y
    vertex_attr(g, "room") <- "basic"
    vertex_attr(g, "room", left_entry) <- "entry_left"
    vertex_attr(g, "room", right_entry) <- "entry_right"
    
    g <- delete_vertices(g, drop)
    return(g)
}


.carve_corners <- function(size, y, x, h, w) {
    bl <- do.call(
        rbind,
        lapply(seq_len(size), function(i) data.frame(
            y = round(y[i]),
            x = seq_len(floor(x[i]) + 1) - 1
        ))
    )
    br <- do.call(
        rbind,
        lapply(seq_len(size), function(i) data.frame(
            y = round(y[i]),
            x = w - seq_len(floor(x[i]) + 1)
        ))
    )
    tl <- do.call(
        rbind,
        lapply(seq_len(size), function(i) data.frame(
            y = h - ceiling(y[i]),
            x = seq_len(floor(x[i]) + 2) - 1
        ))
    )
    tr <- do.call(
        rbind,
        lapply(seq_len(size), function(i) data.frame(
            y = h - ceiling(y[i]),
            x = w + 1 - seq_len(ceiling(x[i]) + 2)
        ))
    )
    lr <- expand.grid(y = seq_len(h) - 1L, x = c(0, w - 1L))
    
    corners <- unique(rbind(bl, br, tl, tr, lr))
    return(corners)
}

.carve_rooms <- function(g, n = 5) {
    for(i in seq_len(n) ) {
        g <- .carve_room(g)
    }
    return(g)
}

.carve_room <- function(g) {
    vy <- vertex_attr(g, "y")
    vx <- vertex_attr(g, "x")
    nonbasic <- vertex_attr(g, "room") != "basic"
    lower_empty <- any(nonbasic & (vy < max(vy) / 2))
    lower_adj <- as.numeric(lower_empty & vy < max(vy) / 2)[!nonbasic]
    
    gg <- delete_vertices(g, nonbasic)
    nnn <- neighborhood_size(gg, 3L)
    if( any(nnn == 25L) & !lower_empty ){
        nn <- nnn
    }else{
        nn <- neighborhood_size(gg, 2L) + lower_adj
    }
    if(!any(nn %in% c(13L, 25L))) return(g)
    
    target <- sample(which(nn == max(nn)), 1)
    
    ty <- vertex_attr(gg, "y")[target]
    tx <- vertex_attr(gg, "x")[target]
    target <- which(vy == ty & vx == tx)
    # determine which direction to grow rest of room; move away from the center
    x_adj <- if( tx < mean( range(vx) ) ) -1L else 1L
    y_adj <- if( ty < mean( range(vy) * 2 / 3 ) ) -1L else 1L
    # except if you're at the border.
    if(tx %in% range(vx)) x_adj <- -x_adj
    if(ty %in% range(vy)) y_adj <- -y_adj
    
    drop <- vy == ty & vx == tx + x_adj |
        vy == ty + y_adj & vx == tx  |
        vy == ty + y_adj & vx == tx + x_adj
    vertex_attr(g, "room")[target] <- paste0(
        if(y_adj > 0) "b" else "t", if(x_adj > 0) "l" else "r"
    )
    g <- delete_vertices(g, drop)
    return(g)
}

.randomize_weights <- function(g) set_edge_attr(
    g, "weight", value = exp(rnorm(length(E(g))))
)

carve_maze <- function(g) mst(.randomize_weights(g))

hex_lines <- function(size) {
    size <- floor(size)
    # height, width, radius, side (apothem)
    h <- size * 2
    w <- round(h * sqrt(3) / 2)
    # zero inclusive, origin at 0,0
    w <- w - 1
    h <- h - 1
    
    r <- h / 2
    s <- w / 2
    
    data.frame(
        x = c(s, w, w, s, 0, 0, s),
        y = c(0, r / 2, h - r / 2, h, h - r / 2, r / 2,  0)
    )
}

# Define plot params ----

size <- 8
colour_path <- "#FEC44F"
colour_wall <- "#EC7014"
colour_bg <- "#FE9929"
colour_border <- "darkred"
colour_text <- "dodgerblue"

# https://fontlibrary.org/en/font/dancing-font
register_font(
    name = "Dancing Script",
    plain = paste0(asset_dir, "DancingScript-Regular.ttf"),
    bold = paste0(asset_dir, "DancingScript-Bold.ttf")
)

# https://openmoji.org/
register_font(
    name = "OpenMoji",
    plain = paste0(asset_dir, "OpenMoji-color-colr0_svg.ttf")
)

label_content <- c(
    "🧬", "🍒", "🍄", "🐉", "🕷", "🫀", "🧠", "🫁", "🦴", "🦠", "💩", "💊",
    "👽", "🍩", "💉", "🦟", "🧪"
)

plot_graph <- hex_graph(size) |>
    .carve_rooms() |>
    carve_maze() |>
    as_tbl_graph()

room <- plot_graph |> pull(room)
label <- rep("", length(room))

label[room %in% c("bl", "tl", "br", "tr")] <- sample(
    label_content, sum(room %in% c("bl", "tl", "br", "tr"))
)

plot_graph <- plot_graph |> mutate(label = label)

# ggplot2 code ----

hex_plot <- plot_graph |>
    ggraph(layout = "manual", y = y, x = x) +
    geom_polygon(data = hex_lines(size), aes(x = x, y = y), fill = colour_bg) +
    geom_edge_link(linewidth = 34 / 6, colour = colour_wall) +
    geom_node_tile(
        aes(filter = room %in% c("basic")),
        linewidth = 2, width = 2 / 3, height = 2 / 3,
        fill = colour_path, colour = colour_wall,
        linejoin = "round", lineend = "round"
    ) +
    geom_node_tile(
        aes(x = x - 0.5 + grepl("l$", room), y = y +0.5 - grepl("^t", room),
            filter = room %in% c("bl", "br", "tl", "tr")),
        linewidth = 2, colour = colour_wall, fill = colour_path, width = 11 / 6,
        height = 11 / 6, linejoin = "round", lineend = "round"
    ) +
    geom_node_point(
        aes(alpha = room %in% c("entry_left", "entry_right")),
        size = 5, shape = 22, fill = colour_path, colour = colour_wall,
        stroke = 11 / 6
    ) +
    # Title box
    annotate(
        "tile",
        x = round(size * sqrt(3) - 1) / 2,
        y = mean(round(c(0, 1) + size * 2 / 3)),
        height = 10 / 6, width = 35 / 6, linewidth = 2,
        fill = colour_path, colour = colour_wall,
        linejoin = "round", lineend = "round"
    ) +
    geom_edge_link(linewidth = 2, colour = colour_path) +
    geom_path(
        data = hex_lines(size), aes(x = x, y = y),
        linewidth = 6, linejoin = "round", lineend = "round",
        colour = colour_border
    ) +
    geom_node_text(
        aes(x = x - 0.5 + grepl("l$", room), y = y + 0.5 - grepl("^t", room),
            filter = room %in% c("bl", "br", "tl", "tr"), label = label
        ), family = "OpenMoji", size = 16) +
    annotate(
        "text",
        label = "ariadne",
        x = round(size * sqrt(3) - 1) / 2,
        y = mean(round(c(0, 1) + size * 2 / 3)),
        hjust = 1 / 2, vjust = 1 / 2,
        size = 60, size.unit = "pt",
        family = "Dancing Script",
        fontface = "bold",
        colour = colour_border
    ) +
    scale_alpha_identity() +
    coord_equal() +
    theme_void()

# Get sticker! ----

file_path <- paste0(asset_dir, "ariadne_logo.png")
ggsave(file_path, width = 10, height = 10, units = "in", scale = 1)
