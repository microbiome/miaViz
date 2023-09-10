test_that("plot RDA/CCA", {
  data("Tengeler2020", package = "mia")
  tse <- Tengeler2020
  
  ### 1). TEST error messages ###
  
  # Object without reducedDim
  expect_error(plotRDA(tse), 'argument "dimred" is missing, with no default')
  expect_error(plotRDA(tse, "RDA"), "'dimred' must specify reducedDim.")
  
  # Run/calculate RDA
  tse <- runRDA(tse, assay ~ patient_status + cohort)
  rda <- calculateRDA(tse, assay ~ patient_status + cohort)
  
  # Minimal functionality
  expect_no_error(plotRDA(tse, "RDA"))

  # Wrong-entry scenarios
  expect_error(plotRDA(tse, "RDA", colour_by = "wrong colname"),
               "'colour_by' must match the name of a column in colData.")
  expect_error(plotRDA(tse, "RDA", colour_by = "cohort", shape_by = "wrong colname"),
               "'shape_by' must match the name of a column in colData.")
  expect_error(plotRDA(tse, "RDA", add.ellipse = "invalid entry"),
               "'add.ellipse' must be one of c(TRUE, FALSE, 'fill', 'color', 'colour').",
               fixed = TRUE)
  ## add more tests here

  ### 2). TEST plot layers ###
  
  el_true <- plotRDA(tse, "RDA", colour_by = "patient_status")
  el_false <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = FALSE)
  el_col <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = "colour")
  el_fill <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = "fill")
  
  # filled ellipse has one more layer than no ellipse plot
  expect_equal(length(ggplot_build(el_true)[["data"]]), 4)
  expect_equal(length(ggplot_build(el_false)[["data"]]), 3)
  
  # coloured ellipse but not filled ellipse plot has all 0 alpha values
  expect_true(all(ggplot_build(el_col)[["data"]][[2]][["alpha"]] == 0))
  expect_false(all(ggplot_build(el_fill)[["data"]][[2]][["alpha"]] == 0))
  
  p_aes <- plotRDA(tse, "RDA", colour_by = "patient_status", ellipse.alpha = 0.5,
                   ellipse.linewidth = 0.2, ellipse.linetype = 3)
  
  # ellipse aesthetics are correctly defined in ggplot
  expect_true(all(ggplot_build(p_aes)[["data"]][[2]][["alpha"]] == 0.5))
  expect_true(all(ggplot_build(p_aes)[["data"]][[2]][["linewidth"]] == 0.2))
  expect_true(all(ggplot_build(p_aes)[["data"]][[2]][["linetype"]] == 3))
  
  ## add more tests here
  
})
