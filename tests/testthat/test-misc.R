context("DMN")
test_that("DMN", {
    data(dmn_se)
    # plot the fit
    actual <- plotDMNFit(dmn_se, type = "laplace")
    expect_s3_class(actual, c("gg","ggplot"))
})
