
context("plot tree")
test_that("plot tree", {
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .get_object_and_trimmed_tree
    expect_error(miaViz:::.get_object_and_trimmed_tree())
    actual <- miaViz:::.get_object_and_trimmed_tree(x["549322",])
    expect_s3_class(rowTree(actual),"phylo")
    expect_s4_class(actual,"TreeSummarizedExperiment")
    expect_equal(unique(rowTree(actual)$tip.label), c("549322"))
    actual <- miaViz:::.get_object_and_trimmed_tree(x)
    expect_equal(rowTree(actual)$tip.label, rownames(x))
    actual <- miaViz:::.get_object_and_trimmed_tree(x, relabel = TRUE)
    expect_equal(rowTree(actual)$tip.label[1L], "Class:Thermoprotei")
    #
    library(scater)
    library(mia)
    data(GlobalPatterns)
    altExp(GlobalPatterns,"genus") <- agglomerateByRank(GlobalPatterns,"Genus", make_unique = FALSE)
    altExp(GlobalPatterns,"genus") <- addPerFeatureQC(altExp(GlobalPatterns,"genus"))
    rowData(altExp(GlobalPatterns,"genus"))$log_mean <- log(rowData(altExp(GlobalPatterns,"genus"))$mean)
    top_taxa <- getTop(altExp(GlobalPatterns,"genus"),
                           method="mean",
                           top=100L,
                           assay.type="counts")
    #
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected"))
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    # plot with tip labels
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       show.label = TRUE))
    expect_true(all(c("colour_by") %in% colnames(plot$data)))
    # plot with selected labels
    labels <- c("Genus:Providencia", "Genus:Morganella", "0.961.60")
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected",
                                       show.label = labels,
                                       layout="rectangular"))
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    # Test that error occurs if tree.name is wrong
    expect_error( plotRowTree(GlobalPatterns, tree.name = "test") )
    expect_error( plotRowTree(GlobalPatterns, tree.name = NULL) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = c("test", "phylo")) )
    expect_error( plotColTree(GlobalPatterns, tree.name = 1) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = TRUE) )
})

# Helpers for heatmap tests
.hm_test_data <- function(n = 20L){
    data(GlobalPatterns, package = "mia")
    x <- mia::agglomerateByRank(GlobalPatterns, "Genus")
    x <- x[mia::getTop(x, top = n), ]
    set.seed(1)
    rowData(x)$a <- rnorm(nrow(x))
    rowData(x)$b <- rnorm(nrow(x))
    rowData(x)$pos <- abs(rowData(x)$a)
    rowData(x)$cat <- rep(c("yes", "no"), length.out = nrow(x))
    rowData(x)$fct <- factor(
        rep(c("low", "mid", "high"), length.out = nrow(x)),
        levels = c("low", "mid", "high"))
    rowData(x)$g1 <- rep(c(TRUE, FALSE), length.out = nrow(x))
    rowData(x)$g2 <- rep(c(FALSE, TRUE, NA), length.out = nrow(x))
    x
}
.hm_layer <- function(p, geom = "GeomTile"){
    which(vapply(
        p$layers, function(l) class(l$geom)[1L] == geom, logical(1L)))
}
.hm_tiles <- function(p) p$layers[[.hm_layer(p)]]$data
.hm_fills <- function(p){
    ggplot2::ggplot_build(p)$data[[.hm_layer(p)]]
}
.hm_limits <- function(p){
    ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")$get_limits()
}

test_that("plotRowTree heatmap: numeric values are aligned with tips", {
    x <- .hm_test_data()
    # No heatmap, no extra layers
    p <- plotRowTree(x)
    expect_length(.hm_layer(p), 0L)
    expect_length(.hm_layer(p, "GeomSegment"), 0L)
    for( layout in c("circular", "rectangular", "fan", "dendrogram") ){
        for( order.tree in c(FALSE, TRUE) ){
            p <- suppressWarnings(plotRowTree(
                x, heatmap.by = c("a", "b"), layout = layout,
                order.tree = order.tree, show.label = rownames(x)))
            expect_s3_class(p, "ggtree")
            expect_length(.hm_layer(p), 1L)
            hm <- .hm_tiles(p)
            tips <- p$data[p$data$isTip, ]
            expect_equal(nrow(hm), 2L * nrow(x))
            expect_equal(levels(hm$variable), c("a", "b"))
            # Each tile is on the same row as its tip
            for( var in c("a", "b") ){
                sub <- hm[hm$variable == var, ]
                expect_equal(sub$y, tips$y)
                expect_equal(
                    sub$value, rowData(x)[[var]][match(tips$label, rownames(x))])
            }
            # Columns are outside the tree, in the given order
            expect_true(all(hm$x > max(tips$x)))
            expect_true(
                unique(hm$x[hm$variable == "a"]) <
                unique(hm$x[hm$variable == "b"]))
            # Tip labels are outside the heatmap
            lab <- ggplot2::ggplot_build(p)$data
            lab <- lab[[length(lab)]]
            expect_true(all(lab$x > max(hm$x)))
        }
    }
})

test_that("plotRowTree heatmap: numeric colour scales", {
    x <- .hm_test_data()
    # Negative values -> diverging scale where zero is the mid colour
    rowData(x)$z <- rowData(x)$a
    rowData(x)$z[1:3] <- 0
    p <- plotRowTree(x, heatmap.by = "z")
    fills <- .hm_fills(p)
    zero <- .hm_tiles(p)$value == 0
    expect_true(all(fills$fill[zero] == "#FFFFFF"))
    sc <- p$scales$get_scales("fill")
    expect_equal(sc$limits, c(-1, 1) * max(abs(rowData(x)$z)))
    # Custom diverging colours: extreme value gets the high/low colour
    p <- plotRowTree(
        x, heatmap.by = "z", heatmap.colours = c("blue", "white", "red"))
    fills <- .hm_fills(p)
    vals <- .hm_tiles(p)$value
    m <- max(abs(vals))
    if( any(vals == m) ) expect_equal(fills$fill[vals == m][1], "#FF0000")
    if( any(vals == -m) ) expect_equal(fills$fill[vals == -m][1], "#0000FF")
    # Positive values only -> sequential scale by default, diverging on request
    p <- plotRowTree(x, heatmap.by = "pos")
    expect_null(p$scales$get_scales("fill")$limits)
    p <- plotRowTree(x, heatmap.by = "pos", heatmap.diverging = TRUE)
    expect_equal(
        p$scales$get_scales("fill")$limits, c(-1, 1) * max(rowData(x)$pos))
    # Negative values with diverging = FALSE -> sequential scale
    p <- plotRowTree(x, heatmap.by = "a", heatmap.diverging = FALSE)
    expect_null(p$scales$get_scales("fill")$limits)
    # Custom sequential colours
    p <- plotRowTree(x, heatmap.by = "pos", heatmap.colours = c("white", "red"))
    fills <- .hm_fills(p)
    vals <- .hm_tiles(p)$value
    expect_equal(fills$fill[which.max(vals)], "#FF0000")
    expect_equal(fills$fill[which.min(vals)], "#FFFFFF")
    # Wrong number of colours
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.colours = "red"))
    expect_error(plotRowTree(
        x, heatmap.by = "pos", heatmap.colours = c("a", "b", "c")))
    # Values from assay
    p <- plotRowTree(x, heatmap.by = colnames(x)[1:2], assay.type = "counts")
    hm <- .hm_tiles(p)
    tips <- p$data[p$data$isTip, ]
    expect_equal(
        hm$value[hm$variable == colnames(x)[1]],
        unname(assay(x, "counts")[match(tips$label, rownames(x)), 1]))
})

test_that("plotRowTree heatmap: categorical and logical values", {
    x <- .hm_test_data()
    # Factor levels are kept
    p <- plotRowTree(x, heatmap.by = "fct")
    expect_equal(levels(.hm_tiles(p)$value), c("low", "mid", "high"))
    expect_equal(.hm_limits(p), c("low", "mid", "high"))
    # Character values from several columns share one scale
    rowData(x)$cat2 <- rep(c("maybe", "no"), length.out = nrow(x))
    p <- plotRowTree(x, heatmap.by = c("cat", "cat2"))
    expect_setequal(levels(.hm_tiles(p)$value), c("maybe", "no", "yes"))
    # Named colours are used
    cols <- c(yes = "red", no = "blue")
    p <- plotRowTree(x, heatmap.by = "cat", heatmap.colours = cols)
    fills <- .hm_fills(p)
    vals <- as.character(.hm_tiles(p)$value)
    expect_equal(fills$fill, unname(cols[vals]))
    expect_error(plotRowTree(
        x, heatmap.by = "cat", heatmap.colours = c(yes = "red")))
    expect_error(plotRowTree(x, heatmap.by = "fct", heatmap.colours = "red"))
    # Logical: TRUE gets the variable name, FALSE and NA are left blank
    p <- plotRowTree(x, heatmap.by = c("g1", "g2"))
    hm <- .hm_tiles(p)
    tips <- p$data[p$data$isTip, ]
    expect_equal(levels(hm$value), c("g1", "g2"))
    ord <- match(tips$label, rownames(x))
    expect_equal(
        !is.na(hm$value[hm$variable == "g1"]), rowData(x)$g1[ord])
    expect_equal(
        !is.na(hm$value[hm$variable == "g2"]), rowData(x)$g2[ord] %in% TRUE)
    expect_equal(.hm_limits(p), c("g1", "g2"))
    p <- plotRowTree(
        x, heatmap.by = c("g1", "g2"),
        heatmap.colours = c(g1 = "red", g2 = "blue"))
    fills <- .hm_fills(p)
    # Absent values are not coloured
    expect_setequal(unique(fills$fill), c("red", "blue", NA))
    # Diverging scale is only for numeric values
    expect_error(plotRowTree(x, heatmap.by = "cat", heatmap.diverging = TRUE))
    # Heatmap fill scale is separate from the tip point fill scale
    p <- plotRowTree(
        x, tip.colour.by = "Phylum", heatmap.by = "cat",
        heatmap.colours = cols)
    b <- ggplot2::ggplot_build(p)
    pts <- b$data[[.hm_layer(p, "GeomPoint")[1]]]
    expect_false(any(pts$fill %in% c("red", "blue")))
    expect_setequal(unique(b$data[[.hm_layer(p)]]$fill), c("red", "blue"))
})

test_that("plotRowTree heatmap: visual options", {
    x <- .hm_test_data()
    # Defaults: guides, white borders and column names
    p <- plotRowTree(x, heatmap.by = c("a", "b"))
    expect_length(.hm_layer(p, "GeomSegment"), 1L)
    expect_equal(p$layers[[.hm_layer(p)]]$aes_params$colour, "white")
    txt <- p$layers[[.hm_layer(p, "GeomText")]]$data
    expect_equal(txt$label, c("a", "b"))
    # Options switched off
    p <- plotRowTree(
        x, heatmap.by = c("a", "b"), heatmap.guides = FALSE,
        heatmap.border = NA, heatmap.colnames = FALSE)
    expect_length(.hm_layer(p, "GeomSegment"), 0L)
    expect_true(is.na(p$layers[[.hm_layer(p)]]$aes_params$colour))
    expect_length(.hm_layer(p, "GeomText"), 0L)
    # Width and offset are relative to tree depth
    p <- plotRowTree(
        x, heatmap.by = c("a", "b"), heatmap.width = 0.5,
        heatmap.offset = 0.1)
    hm <- .hm_tiles(p)
    depth <- diff(range(p$data$x))
    col_width <- 0.5 * depth / 2
    expect_equal(
        min(hm$x), max(p$data$x[p$data$isTip]) + 0.1 * depth + col_width / 2)
    expect_equal(diff(unique(hm$x)), col_width)
    # Legend title
    p <- plotRowTree(x, heatmap.by = "a")
    expect_equal(p$scales$get_scales("fill")$name, "a")
    p <- plotRowTree(x, heatmap.by = c("a", "b"))
    expect_equal(p$scales$get_scales("fill")$name, "value")
    p <- plotRowTree(x, heatmap.by = c("a", "b"), heatmap.title = "LFC")
    expect_equal(p$scales$get_scales("fill")$name, "LFC")
    # US spelling
    p <- plotRowTree(x, heatmap.by = "cat",
        heatmap.colors = c(yes = "red", no = "blue"))
    expect_setequal(unique(.hm_fills(p)$fill), c("red", "blue"))
})

test_that("plotColTree heatmap", {
    x <- .hm_test_data(10L)
    colTree(x) <- ape::rtree(ncol(x), tip.label = colnames(x))
    p <- suppressWarnings(plotColTree(x, heatmap.by = "SampleType"))
    hm <- .hm_tiles(p)
    tips <- p$data[p$data$isTip, ]
    expect_equal(
        as.character(hm$value),
        as.character(x$SampleType[match(tips$label, colnames(x))]))
})

test_that("plotRowTree heatmap: input checks", {
    x <- .hm_test_data(10L)
    expect_error(plotRowTree(x, heatmap.by = c("a", "cat")))
    expect_error(plotRowTree(x, heatmap.by = c("g1", "cat")))
    expect_error(plotRowTree(x, heatmap.by = c("g1", "a")))
    expect_error(plotRowTree(x, heatmap.by = "unknown"))
    expect_error(plotRowTree(x, heatmap.by = c("a", "a")))
    expect_error(plotRowTree(x, heatmap.by = character(0)))
    expect_error(plotRowTree(x, heatmap.by = 1))
    expect_error(plotRowTree(x, heatmap.by = NA_character_))
    for( layout in c("equal_angle", "daylight", "inward_circular") ){
        expect_error(plotRowTree(x, heatmap.by = "a", layout = layout),
            "not supported")
    }
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.width = -1))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.width = "a"))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.offset = -1))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.colnames = "yes"))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.colnames.angle = "a"))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.font.size = 0))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.title = 1))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.title = c("a", "b")))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.colours = 1:3))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.diverging = "yes"))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.border = c("a", "b")))
    expect_error(plotRowTree(x, heatmap.by = "a", heatmap.guides = NA))
})
