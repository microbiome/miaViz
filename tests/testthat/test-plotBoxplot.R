test_that("plotBoxplot check input tests", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rowData(tse)[["value"]] <- rnorm(nrow(tse))
    colData(tse)[["value"]] <- rnorm(ncol(tse))
    # Error if no variable or assay specified
    expect_error(
        plotBoxplot(tse),
        "Please specify either 'assay.type', 'row.var', or 'col.var'.")
    # Error if invalid feature selection
    expect_error(plotBoxplot(tse, row.var = "value", features = rownames(tse)))
    expect_error(plotBoxplot(tse, col.var = "value", features = rownames(tse)))
    # Error if invalid assay.type
    expect_error(plotBoxplot(tse, assay.type = "test"))
    # Error if invalid row.var
    expect_error(plotBoxplot(tse, row.var = "test"))
    # Error if invalid col.var
    expect_error(plotBoxplot(tse, col.var = "test"))
    # Error if both col.var and assay.type are provided
    expect_error(plotBoxplot(tse, col.var = "value", assay.type = "counts"))
    # Error if both row.var and assay.type are provided
    expect_error(plotBoxplot(tse, row.var = "value", assay.type = "counts"))
    # Error if non-existent feature is selected
    expect_error(plotBoxplot(tse, features = "test", assay.type = "counts"))
    # Error if invalid layout
    expect_error(plotBoxplot(tse, add.box = "test", assay.type = "counts"))
    # Error if facet.by column missing
    expect_error(
        plotBoxplot(tse, assay.type = "counts", facet.by = "missing_batch")
    )
    # Error if fill.by column missing
    expect_error(
        plotBoxplot(tse, assay.type = "counts", fill.by = "missing_batch")
    )
})

test_that("plotBoxplot returns a ggplot object for valid inputs", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rowData(tse)[["value"]] <- rnorm(nrow(tse))
    colData(tse)[["value"]] <- rnorm(ncol(tse))

    # Check for ggplot object when assay.type is specified
    p <- plotBoxplot(tse, assay.type = "counts")
    expect_s3_class(p, "ggplot")

    # Check for ggplot object when row.var is specified
    p <- plotBoxplot(tse, row.var = "value")
    expect_s3_class(p, "ggplot")

    # Check for ggplot object when col.var is specified
    p <- plotBoxplot(tse, col.var = "value", fill.by = "group")
    expect_s3_class(p, "ggplot")

    # Check that the values are correct
    expect_equal(tse[["value"]], p$data$value)
    expect_equal(tse[["group"]], p$data$group)
})

test_that("plotBoxplot correctly handles features argument", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rownames(tse) <- paste0("gene", seq_len(nrow(tse)))

    # Valid feature subsetting
    p <- plotBoxplot(tse, assay.type = "counts", features = "gene1")
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$feature == "gene1"))
})

test_that("plotBoxplot facet.by splits plots correctly", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    colData(tse)$batch <- sample(
        c("batch1", "batch2"), ncol(tse), replace = TRUE)

    p <- plotBoxplot(tse, assay.type = "counts", facet.by = "batch")
    expect_s3_class(p, "ggplot")
    expect_true("batch" %in% colnames(p$data))
})

test_that("plotBoxplot handles NA values in data", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    assay(tse, "counts")[1, 1] <- NA  # introduce NA

    p <- plotBoxplot(tse, assay.type = "counts", point.offset = "jitter")
    expect_s3_class(p, "ggplot")
    expect_true(any(is.na(p$data$counts)))
})

test_that("plotBoxplot can color, fill and specify x based on same column", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    colData(tse)[["value"]] <- sample(c("A", "B"), ncol(tse), replace = TRUE)

    p <- plotBoxplot(
        tse, assay.type = "counts", x = "value", colour.by = "value",
        fill.by = "value", point.offset = "none")
    p_data <- ggplot_build(p)

    # Colour and fill should be the same
    coloring <- p_data@data[[2]][ , c("colour", "fill")]
    expect_equal(coloring[[1]], coloring[[2]])
    # There is 2 groups, so there should be 2 unique rows (because jitter was
    # not applied)
    expect_true(
        p_data@data[[2]][, c("x", "colour", "fill")] |> unique() |> nrow() == 2L
    )
})

test_that("plotBoxplot draws paired lines in time order", {
    tse <- makeTSE(nrow = 1L, ncol = 6L)
    assayNames(tse) <- "counts"
    colData(tse)[["subject"]] <- rep(c("a", "b"), each = 3L)
    colData(tse)[["time"]] <- rep(c("t1", "t2", "t3"), 2L)
    colData(tse)[["score"]] <- c(1, 5, 2, 3, 3, 9)

    p <- plotBoxplot(
        tse, col.var = "score", x = "time", pair.by = "subject",
        add.change = TRUE, point.offset = "none")
    is_path <- vapply(
        p$layers, function(l) inherits(l$geom, "GeomPath"), logical(1L))
    lines <- layer_data(p, which(is_path))

    # Each line follows the time points in order
    for( g in unique(lines[["group"]]) ){
        expect_equal(as.numeric(lines[lines[["group"]] == g, "x"]), c(1, 2, 3))
    }
    # Each segment is colored by the change to the next time point
    expect_equal(
        p$data[["difference"]][p$data[["subject"]] == "a"], c(4, -3, NA))
    expect_equal(
        p$data[["difference"]][p$data[["subject"]] == "b"], c(0, 6, NA))
})

test_that("plotBoxplot: size.by sets the title of the size legend", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    colData(tse)[["value"]] <- rnorm(ncol(tse))
    colData(tse)[["depth"]] <- seq_len(ncol(tse))
    p <- plotBoxplot(tse, col.var = "value", size.by = "depth")
    expect_identical(p$labels$size, "depth")
    expect_null(p$labels$shape)
})

test_that("plotBoxplot works with variables named like internal objects", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    tse[["group"]] <- rep(c("a", "b"), each = 5L)
    tse[["col"]] <- rep(c("c", "d"), 5L)
    vals <- as.numeric(seq_len(ncol(tse)))
    # Columns named "x", "y" or "df" masked internal variables in aes() (#239)
    for( nm in c("y", "x", "df") ){
        tse[[nm]] <- vals
        p <- plotBoxplot(
            tse, col.var = nm, x = "group", colour.by = "col",
            point.offset = "none")
        pts <- layer_data(p, 2L)
        expect_equal(sort(pts[["y"]]), vals)
        expect_identical(length(unique(pts[["colour"]])), 2L)
        tse[[nm]] <- NULL
    }
    # Variables that would be overwritten by internal columns give an error
    tse[["x_point"]] <- vals
    expect_error(
        plotBoxplot(tse, col.var = "x_point", x = "group"), "reserved")
})

test_that("plotBoxplot places significance brackets above all spanned groups", {
    skip_if_not_installed("ggpubr")
    tse <- makeTSE(nrow = 1L, ncol = 12L)
    assayNames(tse) <- "counts"
    colData(tse)[["grp"]] <- rep(c("A", "B", "C"), each = 4L)
    colData(tse)[["batch"]] <- rep(c("b1", "b2"), 6L)
    # The middle group is the highest
    colData(tse)[["score"]] <- c(1:4, 11:14, 1:4)
    pvals <- data.frame(
        group1 = c("A", "A", "B"), group2 = c("B", "C", "C"),
        p = c(0.01, 0.5, 0.01))

    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none")
    pos <- attributes(p$data)[["p.value"]]
    ac <- pos[["group1"]] == "A" & pos[["group2"]] == "C"
    # A-C bracket must be above B, even though A and C are low
    expect_true(pos[["y.position"]][ac] > max(tse[["score"]]))
    # All brackets overlap in x axis, so they must not be at the same height
    expect_equal(length(unique(pos[["y.position"]])), 3L)
    # The widest bracket is the top one
    expect_equal(pos[["y.position"]][ac], max(pos[["y.position"]]))

    # step.increase is applied only once; ggpubr does not shift the brackets
    # further
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none", step.increase = 0.2)
    pos <- attributes(p$data)[["p.value"]]
    is_bracket <- vapply(
        p$layers, function(l) inherits(l$geom, "GeomBracket"), logical(1L))
    br <- layer_data(p, which(is_bracket))
    horizontal <- br[["y"]] == br[["yend"]] & br[["x"]] != br[["xend"]]
    expect_equal(
        sort(unique(br[["y"]][horizontal])), sort(pos[["y.position"]]))
    # The gap between stacked brackets follows step.increase
    y_range <- diff(range(tse[["score"]]))
    expect_equal(diff(sort(pos[["y.position"]])), rep(0.2 * y_range, 2L))

    # With facets, each facet gets its own brackets
    pvals_f <- rbind(
        cbind(pvals, batch = "b1"), cbind(pvals, batch = "b2"))
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals_f,
        facet.by = "batch", point.offset = "none")
    pos <- attributes(p$data)[["p.value"]]
    for( b in c("b1", "b2") ){
        temp <- pos[pos[["batch"]] == b, ]
        ac <- temp[["group1"]] == "A" & temp[["group2"]] == "C"
        b_max <- max(tse[["score"]][tse[["grp"]] == "B" & tse[["batch"]] == b])
        expect_true(temp[["y.position"]][ac] > b_max)
        expect_equal(length(unique(temp[["y.position"]])), 3L)
    }
})

test_that("plotBoxplot adds p-values when groups are given with group.by", {
    skip_if_not_installed("ggpubr")
    tse <- makeTSE(nrow = 1L, ncol = 12L)
    assayNames(tse) <- "counts"
    colData(tse)[["grp"]] <- rep(c("A", "B"), each = 6L)
    colData(tse)[["sex"]] <- rep(c("m", "f"), 6L)
    colData(tse)[["score"]] <- seq_len(12L)
    pvals <- data.frame(
        grp = c("A", "B"), group1 = "f", group2 = "m", p = c(0.1, 0.2))
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", group.by = "sex", p.value = pvals,
        point.offset = "none")
    pos <- attributes(p$data)[["p.value"]]
    expect_true(all(pos[["xmax"]] - pos[["xmin"]] > 0))
    expect_true(all(pos[["y.position"]] > tapply(
        tse[["score"]], tse[["grp"]], max)[pos[["grp"]]]))
})

test_that("plotBoxplot hide.ns drops brackets before placing them", {
    skip_if_not_installed("ggpubr")
    tse <- makeTSE(nrow = 1L, ncol = 12L)
    assayNames(tse) <- "counts"
    colData(tse)[["grp"]] <- rep(c("A", "B", "C"), each = 4L)
    colData(tse)[["score"]] <- c(1:4, 11:14, 1:4)
    # A-C (the widest bracket) is not significant
    pvals <- data.frame(
        group1 = c("A", "A", "B"), group2 = c("B", "C", "C"),
        p = c(0.01, 0.5, 0.01))
    p_all <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none")
    pos_all <- attributes(p_all$data)[["p.value"]]
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none", hide.ns = TRUE)
    pos <- attributes(p$data)[["p.value"]]
    # Only the significant comparisons are placed and drawn
    expect_equal(nrow(pos), 2L)
    expect_false(any(pos[["group1"]] == "A" & pos[["group2"]] == "C"))
    # The hidden bracket leaves no gap: the shown brackets are not higher
    # than the lowest levels used when all brackets are shown
    expect_true(
        max(pos[["y.position"]]) <= sort(pos_all[["y.position"]])[[2L]])
    # Character hide.ns filters on the given column
    pvals[["p.adj"]] <- c(0.01, 0.5, 0.2)
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none", label = "p.adj", hide.ns = "p.adj")
    expect_equal(nrow(attributes(p$data)[["p.value"]]), 1L)
    # Nothing significant: the plot is drawn without brackets
    pvals[["p"]] <- 0.5
    pvals[["p.adj"]] <- NULL
    p <- plotBoxplot(
        tse, col.var = "score", x = "grp", p.value = pvals,
        point.offset = "none", hide.ns = TRUE)
    expect_equal(nrow(attributes(p$data)[["p.value"]]), 0L)
    expect_false(any(vapply(
        p$layers, function(l) inherits(l$geom, c("GeomBracket", "GeomText")),
        logical(1L))))
    expect_no_warning(ggplot2::ggplot_build(p))
    # Invalid values
    expect_error(
        plotBoxplot(tse, col.var = "score", x = "grp", p.value = pvals,
            hide.ns = "q"),
        "hide.ns")
})
