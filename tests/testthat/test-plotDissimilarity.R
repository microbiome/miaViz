test_that("plotDissimilarity", {
    tse <- makeTSE(nrow = 10, ncol = 5)
    colData(tse)[["group"]] <- c("a", "a", "b", "b", "b")
    metadata(tse)[["diss"]] <- dist(t(assay(tse)))
    p <- plotDissimilarity(tse, "diss")
    expect_s3_class(p, "ggplot")
    # One row per pair of samples
    expect_equal(nrow(p$data), 5 * 5)
    # Rows and columns split by sample variable, panels sized by number of
    # samples
    p <- plotDissimilarity(tse, "diss", col.var = "group")
    expect_equal(as.character(p$data$group),
        tse$group[match(p$data$SampleID, colnames(tse))])
    expect_equal(as.character(p$data$group_row),
        tse$group[match(p$data$FeatureID, colnames(tse))])
    expect_equal(p$facet$params$space_free$x, TRUE)
    expect_equal(p$facet$params$space_free$y, TRUE)
    expect_no_error(ggplot2::ggplot_build(p))
    # Without sample names, samples are assumed to be in the order of x
    d <- metadata(tse)[["diss"]]
    attr(d, "Labels") <- NULL
    metadata(tse)[["nolab"]] <- d
    p <- plotDissimilarity(tse, "nolab", col.var = "group")
    expect_setequal(as.character(p$data$SampleID), colnames(tse))
    # Wrong input
    expect_error(plotDissimilarity(tse, "wrong"))
    expect_error(plotDissimilarity(tse, "diss", col.var = "wrong"))
    metadata(tse)[["mat"]] <- as.matrix(metadata(tse)[["diss"]])
    expect_error(plotDissimilarity(tse, "mat"), "dist")
    # Number of samples must match
    expect_error(plotDissimilarity(tse[, 1:4], "diss"), "ncol")
    # Sample names must match
    d <- metadata(tse)[["diss"]]
    attr(d, "Labels") <- paste0("x", 1:5)
    metadata(tse)[["badlab"]] <- d
    expect_error(plotDissimilarity(tse, "badlab"), "colnames")
})

test_that("plotDissimilarity sorts samples", {
    skip_if_not_installed("seriation")
    tse <- makeTSE(nrow = 10, ncol = 6)
    metadata(tse)[["diss"]] <- dist(t(assay(tse)))
    d <- metadata(tse)[["diss"]]
    # Same order on both axes
    p <- plotDissimilarity(tse, "diss")
    ord <- labels(d)[getHeatmapOrder(d, method = "MDS_angle")]
    expect_identical(levels(p$data$SampleID), ord)
    expect_identical(levels(p$data$FeatureID), ord)
    p <- plotDissimilarity(tse, "diss", sort.method = "OLO")
    ord <- labels(d)[getHeatmapOrder(d, method = "OLO")]
    expect_identical(levels(p$data$SampleID), ord)
    # Without sorting, samples follow the order of x
    tse <- tse[, rev(colnames(tse))]
    p <- plotDissimilarity(tse, "diss", sort.method = NULL)
    expect_identical(levels(p$data$SampleID), colnames(tse))
    expect_identical(levels(p$data$FeatureID), colnames(tse))
    expect_error(plotDissimilarity(tse, "diss", sort.method = mean))
})

test_that("plotDissimilarity supports separate row and column sorting", {
    skip_if_not_installed("seriation")
    tse <- makeTSE(nrow = 10, ncol = 6)
    d <- dist(t(assay(tse)))
    metadata(tse)[["diss"]] <- d
    p <- plotDissimilarity(tse, "diss", sort.cols = "OLO")
    expect_identical(
        levels(p$data$FeatureID), labels(d)[getHeatmapOrder(d)])
    expect_identical(
        levels(p$data$SampleID),
        labels(d)[getHeatmapOrder(d, method = "OLO")])
    p <- plotDissimilarity(tse, "diss", sort.rows = NULL)
    expect_identical(levels(p$data$FeatureID), colnames(tse))
    expect_error(plotDissimilarity(tse, "diss", sort.cols = mean), "sort.cols")
})
