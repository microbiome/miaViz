
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
        miaViz:::.get_prevalence_count(1,1,mat, include_lowest = TRUE),49)
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
                                   colour_by = "Phylum")
    expect_s3_class(plot,"ggplot")
    expect_equal(dim(plot$data),c(341,4))
})
