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
  
  # Create a mock dataset
  df <- data.frame(
    Feature = rep(c("Feature1", "Feature2", "Feature3"), times = 2),
    Value = c(2, 4, -1, -3, 5, 7),
    Value_abs = abs(c(2, 4, -1, -3, 5, 7)),
    Sign = c("+", "+", "-", "-", "+", "+"),
    PC = rep(c("PC1", "PC2"), each = 3)
  )
  
  # Create an empty ggplot object for plot_out
  plot_out <- ggplot(df)
  ### 1). TEST: barplot with absolute scale and color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = TRUE, show.color = TRUE)
  expect_s3_class(plot, "ggplot")
  ### 2). TEST: barplot without absolute scale but with color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = FALSE, show.color = TRUE)
  expect_s3_class(plot, "ggplot")
  ### 3). TEST: barplot with absolute scale but no color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = TRUE, show.color = FALSE)
  expect_s3_class(plot, "ggplot")
  ### 4). TEST: lollipop plot with absolute scale and color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "lollipop", absolute.scale = TRUE, show.color = TRUE)
  expect_s3_class(plot, "ggplot")
  ### 5). TEST: lollipop plot without absolute scale but with color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "lollipop", absolute.scale = FALSE, show.color = TRUE)
  expect_s3_class(plot, "ggplot")
  ### 6). TEST: lollipop plot with absolute scale but no color
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "lollipop", absolute.scale = TRUE, show.color = FALSE)
  expect_s3_class(plot, "ggplot")
  ### 7). TEST: error when `absolute.scale` is not a boolean
  expect_error(
    .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = "not boolean", show.color = TRUE),
    "'absolute.scale' must be TRUE or FALSE."
  )
  ### 8). TEST: error when `show.color` is not a boolean
  expect_error(
    .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = TRUE, show.color = "not boolean"),
    "'show.color' must be TRUE or FALSE."
  )
  ### 9). TEST: correct labels in the legend
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = TRUE, show.color = TRUE)
  expect_true("Effect" %in% ggplot_build(plot)$plot$scales$scales[[1]]$name)
  expect_equal(ggplot_build(plot)$plot$scales$scales[[1]]$labels, c("+" = "positive", "-" = "negative"))
  ### 10). TEST: adding sign labels
  plot <- .plot_bar_or_lollipop(plot_out, df, layout = "barplot", absolute.scale = TRUE, show.color = TRUE, show.sign = TRUE)
  expect_s3_class(plot, "ggplot")
})

