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
        add.chance = TRUE, point.offset = "none")
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
