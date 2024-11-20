test_that("plot Eigenvalues", {
    # plotScree handles non-numeric eigenvalues in vector
    expect_error(
        plotScree(c("a", "b", "c")),
        "'x' must be a numeric vector or class 'eigenvals'.")
    # Missing eigenvalues in TreeSummarizedExperiment
    tse <- TreeSummarizedExperiment(
        assays = list(counts = matrix(rpois(1000, 5), ncol = 10)))
    # Add reducedDim without eigenvalues
    reducedDim(tse, "PCA") <- matrix(rnorm(100), ncol = 10)
    expect_error(plotScree(tse, "PCA"))
    # Invalid dimred input
    expect_error(
        plotScree(tse, "invalid_dimred"), 
        "'dimred' must specify a valid reducedDim.")
    
    # Define some eigenvalues for vector-based tests
    eigenvalues <- sort(runif(10), decreasing = TRUE)
    # Check that eigenvalues are plotted from TreeSE or from vector
    attr(reducedDim(tse, "PCA"), "test") <- eigenvalues
    expect_error(plotScree(tse, "PCA"))
    p1 <- plotScree(tse, "PCA", eig.name = "test")
    p2 <- plotScree(eigenvalues)
    # Check if a ggplot object is returned
    expect_s3_class(p1, "ggplot")
    # Check if the plots are equal
    df1 <- ggplot_build(p1)[[1]][[1]]
    df2 <- ggplot_build(p2)[[1]][[1]]
    expect_equal(df1, df2)
    
    # Test with different options
    p <- plotScree(tse, 1, eig.name = "test", add.cumulative = TRUE, n = 10000)
    expect_s3_class(p, "ggplot")
})
