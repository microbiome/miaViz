test_that("plotNMDS", {
    # Test the class of plotNMDS output
    data("esophagus", package = "mia")
    esophagus <- runNMDS(esophagus)
    plot <- plotNMDS(esophagus)
    expect_s3_class(plot,c("gg","ggplot"))
})
