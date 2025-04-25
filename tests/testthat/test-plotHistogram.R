
test_that("plotHistogram check input tssts", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rowData(tse)[["value"]] <- rnorm(nrow(tse))
    colData(tse)[["value"]] <- rnorm(ncol(tse))
    #
    expect_error(
        plotHistogram(tse),
        "Please specify either 'assay.type', 'row.var', or 'col.var'.")
    #
    expect_error(plotHistogram(
        tse, row.var = "value", features = rownames(tse)))
    #
    expect_error(plotHistogram(
        tse, col.var = "value", features = rownames(tse)))
    #
    expect_error(plotHistogram(tse, assay.type = "test"))
    #
    expect_error(plotHistogram(tse, row.var = "test"))
    #
    expect_error(plotHistogram(tse, col.var = "test"))
    #
    expect_error(plotHistogram(tse, col.var = "value", assay.type = "counts"))
    #
    expect_error(plotHistogram(tse, row.var = "value", assay.type = "counts"))
    #
    expect_error(plotHistogram(tse, features = "test", assay.type = "counts"))
    #
    expect_error(plotHistogram(tse, layout = "test", assay.type = "counts"))
    #
    expect_error(plotBarplot(tse, assay.type = "counts"))
})

test_that("plotHistogram returns a ggplot object for valid inputs", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rowData(tse)[["value"]] <- rnorm(nrow(tse))
    colData(tse)[["value"]] <- rnorm(ncol(tse))
    # Check for ggplot object when assay.type is specified
    p <- plotHistogram(tse, assay.type = "counts")
    expect_s3_class(p, "ggplot")

    # Check for ggplot object when row.var is specified
    p <- plotHistogram(tse, row.var = "value")
    expect_s3_class(p, "ggplot")

    # Check for ggplot object when col.var is specified
    p <- plotHistogram(tse, col.var = "value", fill.by = "group")
    expect_s3_class(p, "ggplot")

    # Check that the values are correct
    expect_equal(tse[["value"]], p$data$value)
})

test_that("plotBarplot returns a ggplot object for valid inputs", {
    tse <- makeTSE()
    rowData(tse)[["group"]] <- sample(c("A", "B"), nrow(tse), replace = TRUE)
    # Check for ggplot object when assay.type is specified
    p <- plotBarplot(tse, col.var = "group")
    expect_s3_class(p, "ggplot")
    # Check for ggplot object when row.var is specified
    p <- plotBarplot(tse, row.var = "var1", facet.by = "group")
    expect_s3_class(p, "ggplot")
    # Check that the values are correct
    expect_equal(rowData(tse)[["var1"]], p$data$value)
})
