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
