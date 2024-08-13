test_that("plot Loadings", {
  data("GlobalPatterns", package = "mia")
  tse <- GlobalPatterns
  
  ### 1). TEST error messages ###
  
  # Object without reducedDim
  expect_error(plotLoadings(tse), "No reducedDims found.")
  
  # Run/calculate PCA
  tse <- agglomerateByPrevalence(tse, rank = "Phylum", update.tree = TRUE)
  tse <- transformAssay(tse, method = "clr", pseudocount = 1)
  tse <- scater::runPCA(tse, ncomponents = 5, assay.type = "clr")
  
  # Minimal functionality
  expect_no_error(plotLoadings(tse, dimred = "PCA"))
  
  # Wrong-entry scenarios
  expect_error(plotLoadings(tse, dimred = "PCA", layout = "wrong name"))
  
  # Test that error occurs if tree.name is wrong
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = "test") )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = NULL) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = c("test", "phylo")) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = 1) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = TRUE) )
  
  # Expect that the returned value is always ggplot
  p <- plotLoadings(tse, dimred = "PCA")
  expect_s3_class(p, "ggplot")
  mat <- reducedDim(tse, "PCA")
  p <- plotLoadings(mat, layout = "heatmap")
  expect_s3_class(p, "ggplot")
  p <- plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = TRUE)
  expect_s3_class(p, "ggplot")
})
