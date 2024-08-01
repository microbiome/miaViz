test_that("plot Loadings", {
  data("GlobalPatterns", package = "mia")
  tse <- GlobalPatterns
  
  ### 1). TEST error messages ###
  
  # Object without reducedDim
  expect_error(plotLoadings(tse), "argument \"dimred\" is missing, with no default")
  
  # Run/calculate PCA
  tse <- transformAssay(tse, method = "clr", pseudocount = 1)
  tse <- agglomerateByPrevalence(tse, rank= "Phylum", update.tree = TRUE)
  tse <- scater::runPCA(tse, ncomponents = 5, assay.type = "clr")
  
  # Minimal functionality
  expect_no_error(plotLoadings(tse, dimred = "PCA"))
  
  # Wrong-entry scenarios
  expect_error(plotLoadings(tse, dimred = "PCA", layout = "wrong name"))
  expect_warning(plotLoadings(tse, dimred = "PCA", add.tree= TRUE, layout = "heatmap"),
               "message here")
  
  # Test that error occurs if tree.name is wrong
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = "test") )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = NULL) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = c("test", "phylo")) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = 1) )
  expect_error( plotLoadings(tse, add.tree = TRUE, layout ="heatmap", tree.name = TRUE) )

  
})
