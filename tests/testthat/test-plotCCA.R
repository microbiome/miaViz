test_that("plot RDA/CCA", {
  data("Tengeler2020", package = "mia")
  tse <- Tengeler2020

  ### 1). TEST error messages ###

  # Object without reducedDim
  expect_error(plotRDA(tse), 'argument "dimred" is missing, with no default')
  expect_error(plotRDA(tse, "RDA"), "'dimred' must specify reducedDim.")

  # Run/calculate RDA
  tse <- addRDA(tse, formula = assay ~ patient_status + cohort)
  rda <- getRDA(tse, formula = assay ~ patient_status + cohort)

  # Minimal functionality
  expect_no_error(plotRDA(tse, "RDA"))

  # Wrong-entry scenarios
  expect_error(plotRDA(tse, "RDA", colour_by = "wrong colname"))
  expect_error(plotRDA(tse, "RDA", colour_by = "cohort", shape_by = "wrong colname"))
  expect_error(plotRDA(tse, "RDA", add.ellipse = "invalid value"),
               "'add.ellipse' must be one of c(TRUE, FALSE, 'fill', 'color').",
               fixed = TRUE)
  expect_error(plotRDA(tse, "RDA", add.significance = "invalid value"),
               "'add.significance' must be TRUE or FALSE.")
  expect_error(plotRDA(tse, "RDA", add.expl.var = "invalid value"),
               "'add.expl.var' must be TRUE or FALSE.")
  expect_error(plotRDA(tse, "RDA", repel.labels))
  expect_warning(plotRDA(tse, "RDA", add.significance = FALSE, parse.labels = TRUE),
                 "'parse.labels' was turned off because 'add.significance' is FALSE.")
  expect_warning(plotRDA(tse, "RDA", add.vectors = FALSE),
                 "'add.vectors' is FALSE, so other arguments for vectors and labels will be disregarded.")
  ## add more tests here

  ### 2). TEST plot layers ###

  el_true <- plotRDA(tse, "RDA", colour_by = "patient_status")
  el_false <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = FALSE)
  el_col <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = "colour")
  el_fill <- plotRDA(tse, "RDA", colour_by = "patient_status", add.ellipse = "fill")
  expect_warning(
      vec_false <- plotRDA(tse, "RDA", colour_by = "patient_status", add.vectors = FALSE)
  )
  # Filled ellipse has one more layer than no ellipse plot
  expect_equal(length(ggplot_build(el_true)[["data"]]), 4)
  expect_equal(length(ggplot_build(el_false)[["data"]]), 3)
  # No-vector plot has only 2 layers, (points and ellipse)
  expect_equal(length(ggplot_build(vec_false)[["data"]]), 2)
  # Coloured ellipse but not filled ellipse plot has all 0 alpha values
  expect_true(all(ggplot_build(el_col)[["data"]][[2]][["alpha"]] == 0))
  expect_false(all(ggplot_build(el_fill)[["data"]][[2]][["alpha"]] == 0))

  # Check ggplot aesthetics
  p_aes <- plotRDA(tse, "RDA", colour_by = "patient_status", ellipse.alpha = 0.5,
                   ellipse.linewidth = 0.2, ellipse.linetype = 3, vec.size = 0.6,
                   vec.colour = "red", vec.linetype = 2, arrow.size = 0.15,
                   label.colour = "blue", label.size = 5)
  # Build plot and get data
  p_aes_build <- ggplot_build(p_aes)[["data"]]
  # Ellipse aesthetics are correctly defined in ggplot
  expect_true(all(p_aes_build[[2]][["alpha"]] == 0.5))
  expect_true(all(ggplot_build(p_aes)[[2]][["linewidth"]] == 0.2))
  expect_true(all(ggplot_build(p_aes)[[2]][["linetype"]] == 3))
  # Vector aesthetics are correctly defined in ggplot
  expect_true(all(p_aes_build[[3]][["colour"]] == "red"))
  expect_true(all(p_aes_build[[3]][["linewidth"]] == 0.6))
  expect_true(all(p_aes_build[[3]][["linetype"]] == 2))
  # Label aesthetics are correctly defined in ggplot
  expect_true(all(p_aes_build[[4]][["colour"]] == "blue"))
  expect_true(all(p_aes_build[[4]][["size"]] == 5))
  # Where is arrow size stored?
  # expect_true(arrow_size == 0.15)

  # Vector or label text
  p_vec <- plotRDA(tse, "RDA", colour_by = "patient_status", vec.text = TRUE)
  p_lab <- plotRDA(tse, "RDA", colour_by = "patient_status", vec.text = FALSE)
  # Column "fill" is present in p_vec and missing in p_lab, so length differs by 1
  expect_length(ggplot_build(p_vec)[["data"]][[4]], 29)
  expect_length(ggplot_build(p_lab)[["data"]][[4]], 28)
})
