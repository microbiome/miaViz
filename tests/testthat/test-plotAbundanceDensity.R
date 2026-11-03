test_that("plotAbundanceDensity colours points by a colData variable", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    colData(tse)[["sample_group"]] <- rep(c("a", "b"), length.out = ncol(tse))

    p <- plotAbundanceDensity(tse, colour.by = "sample_group")

    expected <- colData(tse)[p$data[["Sample"]], "sample_group"]
    expect_equal(p$data[["colour_by"]], expected)
})

test_that("plotAbundanceDensity colours points by a rowData variable", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    rowData(tse)[["feature_group"]] <- rep(c("a", "b"), length.out = nrow(tse))

    p <- plotAbundanceDensity(tse, colour.by = "feature_group")

    expected <- rowData(tse)[as.character(p$data[["Y"]]), "feature_group"]
    expect_equal(p$data[["colour_by"]], expected)
})

test_that("plotAbundanceDensity prefers colData over rowData in colour.by", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"
    colData(tse)[["shared"]] <- rep(c("a", "b"), length.out = ncol(tse))
    rowData(tse)[["shared"]] <- rep(c("c", "d"), length.out = nrow(tse))

    p <- plotAbundanceDensity(tse, colour.by = "shared")

    expect_setequal(unique(p$data[["colour_by"]]), c("a", "b"))
})

test_that("plotAbundanceDensity rejects an unknown colour.by", {
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    assayNames(tse) <- "counts"

    expect_error(plotAbundanceDensity(tse, colour.by = "unknown"))
})
