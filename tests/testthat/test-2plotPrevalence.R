
context("plot prevalence")
test_that("plot prevalence", {
    library(DelayedArray)
    #
    data(GlobalPatterns)
    se <- GlobalPatterns
    # .get_prevalence_value
    assay(se,"counts") <- DelayedArray(assay(se,"counts"))
    mat <- assay(se,"counts")
    # .get_prevalence_count
    expect_equal(miaViz:::.get_prevalence_count(1,2,mat),0)
    expect_equal(
        miaViz:::.get_prevalence_count(1,1,mat, include.lowest = TRUE),49)
    # .get_prevalence_plot_data
    actual <- miaViz:::.get_prevalence_plot_data(se,"counts",
                                                 c(0.1,0.2),
                                                 c(0.1,0.2))
    expect_s3_class(actual,"data.frame")
    expect_equal(dim(actual),c(4,3))
    expect_named(actual,c("X","colour_by","Y"))
    #
    plot <- plotPrevalence(se, rank = "Phylum")
    expect_s3_class(plot,"ggplot")
    expect_equal(dim(plot$data),c(70,3))
    plot <- plotPrevalentAbundance(GlobalPatterns, rank = "Family",
                                   colour.by = "Phylum")
    expect_s3_class(plot,"ggplot")
    expect_equal(dim(plot$data),c(341,4))
})

test_that("prevalence plots name the detection thresholds in their unit", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns[seq_len(100), ]
    detection <- c(1, 10, 100)
    p <- plotPrevalence(tse, detection = detection)
    expect_equal(p$labels$x, "Detection threshold (counts)")
    expect_equal(p$labels$y, "Number of features")
    expect_equal(plotPrevalence(tse, xlab = "x")$labels$x, "x")
    # Counts are not shown in percent
    p <- plotRowPrevalence(tse, detection = detection)
    expect_equal(p$labels$x, "Detection threshold (counts)")
    expect_equal(levels(p$data$X), c("1", "10", "100"))
    p <- plotRowPrevalence(tse, detection = NULL)
    expect_equal(p$labels$x, "Detection threshold (relative abundance [%])")
})

test_that("as.relative and as_relative are deprecated", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns[seq_len(100), ]
    tse <- mia::transformAssay(tse, method = "relabundance")
    detection <- c(0.001, 0.01, 0.1)
    # The replacement: relative abundances selected with assay.type
    p <- plotPrevalence(tse, assay.type = "relabundance", detection = detection)
    expect_equal(p$labels$x, "Detection threshold (relabundance)")
    # Both names still work during deprecation and give the same result
    expect_warning(
        p1 <- plotPrevalence(tse, detection = detection, as.relative = TRUE),
        "deprecated")
    expect_warning(
        p2 <- plotPrevalence(tse, detection = detection, as_relative = TRUE),
        "deprecated")
    expect_equal(p1$data, p2$data)
    expect_equal(p1$data$Y, p$data$Y)
    expect_warning(
        plotRowPrevalence(tse, detection = detection, as.relative = TRUE),
        "deprecated")
    expect_warning(
        plotPrevalentAbundance(tse, as.relative = TRUE), "deprecated")
    # No warning without them
    expect_no_warning(plotRowPrevalence(tse, detection = detection))
    expect_no_warning(plotPrevalentAbundance(tse))
})
