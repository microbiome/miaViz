test_that("plot Eigenvalues", {
    data("enterotype", package = "mia")
    tse <- enterotype
    
    tse <- addRDA(
             tse,
             formula = assay ~ ClinicalStatus + Gender + Age,
             FUN = getDissimilarity,
             distance = "bray",
             na.action = na.exclude
             )
    
    # Define some eigenvalues for vector-based tests
    eigenvalues <- sort(runif(10), decreasing = TRUE)
    
    # plotScree handles non-numeric eigenvalues in vector
    expect_error(plotScree(c("a", "b", "c")), 
                 "'x' must be a numeric vector.")
    
    # missing eigenvalues in SingleCellExperiment
    sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(1000, 5), 
                                                              ncol = 10)))
    
    # Add reducedDim without eigenvalues
    reducedDim(sce, "PCA") <- matrix(rnorm(100), ncol = 10)
    
    expect_error(plotScree(sce, "PCA"), 
                 "No eigenvalues found in the specified reducedDim.")
    
    # invalid dimred input in SingleCellExperiment
    expect_error(plotScree(tse, "invalid_dimred"), 
                 "'dimred' must specify a valid reducedDim.")
    
    p <- plotScree(eigenvalues)
    
    # Check if a ggplot object is returned
    expect_s3_class(p, "ggplot")
    
    
    p <- plotScree(tse, "RDA")
    
    # Check if a ggplot object is returned
    expect_s3_class(p, "ggplot")
    })