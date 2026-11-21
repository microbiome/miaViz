
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
    p <- plotPrevalence(tse, detection = detection, as.relative = TRUE)
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
