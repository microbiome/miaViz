.get_boxes <- function(x, y, label, angle, size, hjust, vjust){
    angle <- angle %% 360
    # an approximation?
    width <- nchar(label) * (size * ggplot2:::.pt / 1.4)
    height <- 0.025 * size
    
    hjust <- width * hjust
    vjust <- height * vjust
    
    x1 <- x - hjust
    y1 <- y + vjust
    x2 <- x + hjust
    y2 <- y + vjust
    x3 <- x + hjust
    y3 <- y - vjust
    x4 <- x - hjust
    y4 <- y - vjust
    return(list(x1 = x1,y1 = y1,x2 = x2,y2 = y2,x3 = x3,y3 = y3,x4 = x4,y4 = y4))
}

.test_boxes <- function(boxes){
    box_cmbns <- seq_len(nrow(boxes))
    box_cmbns <- expand.grid(one = box_cmbns, two = box_cmbns)
    box_cmbns <- box_cmbns[box_cmbns$one != box_cmbns$two,]
    edge_cmbns <- data.frame(one = c(1,4,4,3), two = c(2,3,1,2))
    edges <- apply(edge_cmbns, 1L,
                   function(cmbn){
                       .get_edge(boxes[paste0("x",cmbn[1L])],boxes[paste0("y",cmbn[1L])],
                                 boxes[paste0("x",cmbn[2L])],boxes[paste0("y",cmbn[2L])])
                   })
    lapply(edges,
           function(edge){
               
               .sat_test(boxes$x4, boxes$y4, edge)
           })
}

.get_edge <- function(x1, y1, x2, y2){
    f <- x2 == x1
    m <- vector(mode="numeric", length(x1))
    b <- vector(mode="numeric", length(x1))
    m[f] <- Inf
    b[f] <- 0
    m[!f] <- (y2[!f] - y1[!f]) / (x2[!f] - x1[!f])
    b[!f] <- y1[!f] - (m[!f] * x1[!f])
    function(x){
        m * x + b
    }
}

.rotate_x_y <- function(x, y){
    return(data.frame(x = -y, y = x))
}

.sat_test <- function(x, y, FUN){
    y_edge <- FUN(x)
    points_rotated <- .rotate_x_y(x, y_edge)
    side <- points_rotated$x * (x - x) +
        points_rotated$y * (y - y_edge)
    sign(side)
}
