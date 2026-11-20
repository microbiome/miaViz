test_that("plot Loadings", {
  data("GlobalPatterns", package = "mia")
  tse <- GlobalPatterns

  ### 1). TEST error messages ###

  # Object without reducedDim
  expect_error(plotLoadings(tse), "No reducedDims found.")

  # Run/calculate PCA
  tse <- agglomerateByRank(tse, rank = "Phylum", update.tree = TRUE)
  tse <- subsetBy(
    tse, rows = rowMeans(counts > 0) > 0.5, collapse = TRUE,
    update.tree = TRUE)
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
  p <- plotLoadings(tse, dimred = "PCA", layout = "heatmap", add.tree = FALSE)
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


test_that("plotLoadings sorts features in heatmap", {
    data("GlobalPatterns", package = "mia")
    tse <- agglomerateByRank(GlobalPatterns, rank = "Phylum")
    tse <- subsetBy(tse, rows = rowMeans(counts > 0) > 0.5, collapse = TRUE)
    tse <- transformAssay(tse, method = "clr", pseudocount = 1)
    tse <- scater::runPCA(tse, ncomponents = 5, assay.type = "clr")
    loadings <- attr(reducedDim(tse, "PCA"), "rotation")
    # Features are sorted in decreasing order of the component
    p <- plotLoadings(tse, "PCA", layout = "heatmap", sort.method = "PC2")
    expect_identical(
        levels(p$data[["Feature"]]),
        rownames(loadings)[order(loadings[, "PC2"], decreasing = TRUE)])
    # Features are sorted with seriation method
    p <- plotLoadings(tse, "PCA", layout = "heatmap", sort.method = "PCA_angle")
    ser <- seriation::seriate(loadings, method = "PCA_angle", margin = 1L)
    expect_identical(
        levels(p$data[["Feature"]]),
        rownames(loadings)[seriation::get_order(ser)])
    # Clustering gives an order for all features
    p <- plotLoadings(
        tse, "PCA", layout = "heatmap",
        sort.method = bluster::HclustParam())
    expect_setequal(levels(p$data[["Feature"]]), rownames(loadings))
    # Unsupported values give an error
    expect_error(
        plotLoadings(tse, "PCA", layout = "heatmap", sort.method = "test"))
    expect_error(
        plotLoadings(tse, "PCA", layout = "heatmap", sort.method = mean))
})

# Order in which the components are drawn: facets for barplot/lollipop, x axis
# for heatmap
.shown_components <- function(p){
    b <- ggplot2::ggplot_build(p)
    if( inherits(p$facet, "FacetNull") ){
        return(b$layout$panel_params[[1]]$x$get_labels())
    }
    lay <- b$layout$layout
    facet_col <- setdiff(
        colnames(lay), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    as.character(lay[order(lay$PANEL), facet_col[[1]]])
}

test_that("plotLoadings uses readable component names and titles", {
    set.seed(1)
    mat <- matrix(rnorm(60), 20, 3)
    rownames(mat) <- paste0("taxon", seq_len(nrow(mat)))
    # Components without names are labelled "Component 1", ...
    p <- plotLoadings(mat, ncomponents = 2)
    # Component names are returned as character, as before
    expect_type(p$data[["PC"]], "character")
    expect_identical(.shown_components(p), c("Component 1", "Component 2"))
    # The same applies to the default names "V1", "V2", ...
    colnames(mat) <- c("V1", "V2", "V3")
    p <- plotLoadings(mat, ncomponents = 3)
    expect_identical(
        .shown_components(p), c("Component 1", "Component 2", "Component 3"))
    # Other names are kept and components are in their original order
    colnames(mat) <- c("PC10", "PC2", "PC1")
    p <- plotLoadings(mat, ncomponents = 3)
    expect_identical(.shown_components(p), c("PC10", "PC2", "PC1"))
    p <- plotLoadings(mat, layout = "heatmap", ncomponents = 3)
    expect_type(p$data[["PC"]], "character")
    expect_identical(.shown_components(p), c("PC10", "PC2", "PC1"))
    # Sorting by a component still uses the original names
    colnames(mat) <- c("V1", "V2", "V3")
    p <- plotLoadings(mat, layout = "heatmap", ncomponents = 2,
        sort.method = "V2")
    expect_identical(
        levels(p$data[["Feature"]]),
        rownames(mat)[order(mat[, "V2"], decreasing = TRUE)])
    expect_identical(.shown_components(p), c("Component 1", "Component 2"))
    expect_identical(p$labels$fill, "Loading")
    # Axis titles
    p <- plotLoadings(mat, ncomponents = 2)
    expect_identical(p$labels$x, "Loading (absolute value)")
    expect_null(p$labels$y)
    p <- plotLoadings(mat, ncomponents = 2, absolute.scale = FALSE)
    expect_identical(p$labels$x, "Loading")
})

test_that("plotLoadings omits the sign legend when all signs are equal", {
    set.seed(1)
    mat <- matrix(rnorm(40), 20, 2)
    rownames(mat) <- paste0("taxon", seq_len(nrow(mat)))
    # Both signs: legend is shown
    p <- plotLoadings(mat, ncomponents = 2, n = 20)
    expect_false(is.null(ggplot2::get_guide_data(p, "colour")))
    p <- plotLoadings(mat, layout = "barplot", ncomponents = 2, n = 20)
    expect_false(is.null(ggplot2::get_guide_data(p, "fill")))
    # Only positive loadings: legend is omitted
    mat <- abs(mat)
    p <- plotLoadings(mat, ncomponents = 2)
    expect_null(ggplot2::get_guide_data(p, "colour"))
    p <- plotLoadings(mat, layout = "barplot", ncomponents = 2)
    expect_null(ggplot2::get_guide_data(p, "fill"))
    # The user can add the legend back
    p <- plotLoadings(mat, ncomponents = 2) +
        ggplot2::guides(colour = "legend")
    expect_false(is.null(ggplot2::get_guide_data(p, "colour")))
})
