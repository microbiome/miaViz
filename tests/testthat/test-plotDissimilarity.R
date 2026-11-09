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
