test_that("getHeatmapOrder works", {
    skip_if_not_installed("seriation")
    data(GlobalPatterns, package = "mia")
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Phylum")
    tse <- transformAssay(tse, method = "clr", pseudocount = TRUE)
    mat <- assay(tse, "clr")

    # Seriation
    res <- getHeatmapOrder(tse, assay.type = "clr")
    expect_type(res, "integer")
    expect_setequal(res, seq_len(nrow(tse)))
    ref <- seriation::seriate(mat, method = "PCA_angle", margin = 1L) |>
        seriation::get_order() |> unname()
    expect_identical(res, ref)
    # Matrix and SE give the same result
    expect_identical(getHeatmapOrder(mat), res)

    # Columns
    res <- getHeatmapOrder(mat, method = "Heatmap", MARGIN = "cols")
    expect_setequal(res, seq_len(ncol(tse)))
    ref <- seriation::seriate(t(mat), method = "Heatmap", margin = 1L) |>
        seriation::get_order() |> unname()
    expect_identical(res, ref)

    # Single variable
    expect_identical(
        getHeatmapOrder(mat, method = "CL3"),
        order(mat[, "CL3"], decreasing = TRUE))
    expect_identical(
        getHeatmapOrder(mat, method = "CL3", decreasing = FALSE),
        order(mat[, "CL3"]))

    # Clustering
    skip_if_not_installed("bluster")
    res <- getHeatmapOrder(mat, method = bluster::HclustParam())
    expect_setequal(res, seq_len(nrow(tse)))

    # Errors
    expect_error(getHeatmapOrder(mat, method = mean))
    expect_error(getHeatmapOrder(mat, method = c("PCA", "Heatmap")))
    expect_error(getHeatmapOrder(mat, method = "test"))
    expect_error(getHeatmapOrder(mat, MARGIN = 3))
    expect_error(getHeatmapOrder(mat, method = "CL3", decreasing = NULL))
    expect_error(getHeatmapOrder(tse, assay.type = "test"))
})

test_that("getHeatmapOrder works with dist", {
    skip_if_not_installed("seriation")
    mat <- matrix(c(1, 2, 8, 3, 9, 1, 4, 4, 5, 2, 7, 6), nrow = 6,
        dimnames = list(paste0("s", 1:6), NULL))
    d <- dist(mat)
    # Seriation method for distances
    ref <- seriation::seriate(d, method = "MDS_angle") |>
        seriation::get_order() |> unname()
    expect_identical(getHeatmapOrder(d), ref)
    ref <- seriation::seriate(d, method = "OLO") |>
        seriation::get_order() |> unname()
    expect_identical(getHeatmapOrder(d, method = "OLO"), ref)
    # Matrix methods are applied to the distance matrix
    expect_identical(
        getHeatmapOrder(d, method = "PCA_angle"),
        getHeatmapOrder(as.matrix(d), method = "PCA_angle"))
    expect_identical(
        getHeatmapOrder(d, method = "s1"),
        order(as.matrix(d)[, "s1"], decreasing = TRUE))
    expect_error(getHeatmapOrder(d, method = mean))
})
