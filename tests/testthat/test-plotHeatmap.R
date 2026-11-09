test_that("plotHeatmap works for SummarizedExperiment", {

    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    p <- plotHeatmap(
        tse,
        assay.type = "counts"
    )

    expect_s3_class(p, "ggplot")
})

test_that("plotHeatmap supports scaling", {

    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            scale = TRUE
        ),
        "ggplot"
    )

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            center = TRUE
        ),
        "ggplot"
    )

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            scale = TRUE,
            center = TRUE
        ),
        "ggplot"
    )
})

test_that("plotHeatmap supports row and column facets", {

    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            col.var = "SampleType"
        ),
        "ggplot"
    )

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            row.var = "Phylum"
        ),
        "ggplot"
    )

    expect_s3_class(
        plotHeatmap(
            tse,
            assay.type = "counts",
            row.var = "Phylum",
            col.var = "SampleType"
        ),
        "ggplot"
    )
})

test_that("plotHeatmap validates arguments", {

    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "foo"
        )
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            scale = 1
        ),
        "scale"
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            center = 1
        ),
        "center"
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            row.var = "foo"
        )
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            col.var = "foo"
        )
    )
})

test_that("plotHeatmap works for TreeSummarizedExperiment", {

    data(GlobalPatterns)

    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    p <- plotHeatmap(
        tse,
        assay.type = "counts"
    )

    expect_s3_class(p, "ggplot")
})

test_that("plotHeatmap adds row tree", {

    data(GlobalPatterns)

    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    p <- plotHeatmap(
        tse,
        assay.type = "counts",
        show.tree = TRUE
    )

    expect_s3_class(p, "patchwork")
})

test_that("plotHeatmap validates tree arguments", {

    data(GlobalPatterns)

    tse <- GlobalPatterns

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            show.tree = 1
        ),
        "show.tree"
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            show.tree = TRUE,
            row.var = "Phylum"
        ),
        "row.var"
    )

    expect_error(
        plotHeatmap(
            tse,
            assay.type = "counts",
            show.tree = TRUE,
            tree.name = "foo"
        )
    )
})

test_that("Centering centers each feature to zero", {

    data(GlobalPatterns)

    tse <- agglomerateByRank(GlobalPatterns, rank = "Family")

    df <- .get_heatmap_data(
        tse = tse,
        assay.type = "counts",
        center = TRUE,
        scale = FALSE
    )

    means <- tapply(
        df$counts,
        df$FeatureID,
        mean
    )

    expect_true(
        all(abs(means) < 1e-10)
    )
})

test_that("Scaling and centering produce zero mean and unit variance", {

    data(GlobalPatterns)

    tse <- agglomerateByRank(GlobalPatterns, rank = "Family")

    df <- .get_heatmap_data(
        tse = tse,
        assay.type = "counts",
        center = TRUE,
        scale = TRUE
    )

    means <- tapply(
        df$counts,
        df$FeatureID,
        mean
    )

    sds <- tapply(
        df$counts,
        df$FeatureID,
        sd
    )

    expect_true(
        all(abs(means) < 1e-10)
    )

    expect_true(
        all(abs(sds - 1) < 1e-10)
    )
})

test_that("No transformation leaves assay values unchanged", {

    data(GlobalPatterns)

    tse <- agglomerateByRank(GlobalPatterns, rank = "Family")

    df <- .get_heatmap_data(
        tse = tse,
        assay.type = "counts"
    )

    expected <- meltSE(
        tse,
        assay.type = "counts"
    )

    expect_equal(
        df$counts,
        expected$counts,
        check.attributes = FALSE
    )
})

test_that("row tree and heatmap use the same feature order", {

    data(GlobalPatterns)

    tse <- agglomerateByRank(GlobalPatterns, rank = "Family")

    p <- plotHeatmap(
        tse,
        assay.type = "counts",
        show.tree = TRUE,
        show.label = TRUE
    )

    tree_plot <- p[[1]]
    heatmap_plot <- p[[2]]

    # Tip order in the tree
    tree_build <- ggplot2::ggplot_build(tree_plot)
    tree_tips <- tree_build$data[[5]][, c("y", "label")]
    tree_tips <- tree_tips[order(tree_tips$y), ]

    # Feature order in the heatmap
    heatmap_build <- ggplot2::ggplot_build(heatmap_plot)

    heatmap_data <- unique(
        heatmap_build$plot$data[, c("FeatureID")]
    )

    heatmap_order <- as.character(heatmap_data$FeatureID)

    expect_identical(
        heatmap_order,
        tree_tips$label
    )

})

test_that("plotHeatmap sorts rows and columns", {
    skip_if_not_installed("seriation")
    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")

    # By default, the order of the object is kept
    p <- plotHeatmap(tse, assay.type = "counts")
    expect_identical(levels(p$data$FeatureID), rownames(tse))
    expect_identical(levels(p$data$SampleID), colnames(tse))

    # Same method sorts both rows and columns
    p <- plotHeatmap(
        tse, assay.type = "counts", scale = TRUE, center = TRUE,
        sort.method = "PCA_angle")
    mat <- t(scale(t(assay(tse, "counts"))))
    mat[is.nan(mat)] <- 0
    row_order <- getHeatmapOrder(mat, method = "PCA_angle")
    col_order <- getHeatmapOrder(mat, method = "PCA_angle", MARGIN = "cols")
    expect_identical(levels(p$data$FeatureID), rownames(tse)[row_order])
    expect_identical(levels(p$data$SampleID), colnames(tse)[col_order])

    expect_error(plotHeatmap(tse, assay.type = "counts", sort.method = mean))
})

test_that("plotHeatmap sorts only columns when tree is shown", {
    skip_if_not_installed("seriation")
    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")
    p_tree <- plotHeatmap(tse, assay.type = "counts", show.tree = TRUE)
    p <- plotHeatmap(
        tse, assay.type = "counts", show.tree = TRUE, sort.method = "Heatmap")
    # Rows follow the tree
    expect_identical(
        levels(p[[2]]$data$FeatureID), levels(p_tree[[2]]$data$FeatureID))
    col_order <- getHeatmapOrder(
        tse, assay.type = "counts", method = "Heatmap", MARGIN = "cols")
    expect_identical(levels(p[[2]]$data$SampleID), colnames(tse)[col_order])
})

test_that("plotHeatmap supports separate row and column sorting", {
    skip_if_not_installed("seriation")
    data(GlobalPatterns)
    tse <- agglomerateByPrevalence(GlobalPatterns, rank = "Class")
    mat <- as.matrix(assay(tse, "counts"))
    # sort.rows and sort.cols override sort.method
    p <- plotHeatmap(
        tse, assay.type = "counts", sort.method = "PCA_angle",
        sort.cols = "Heatmap")
    row_order <- getHeatmapOrder(mat, method = "PCA_angle")
    col_order <- getHeatmapOrder(mat, method = "Heatmap", MARGIN = "cols")
    expect_identical(levels(p$data$FeatureID), rownames(tse)[row_order])
    expect_identical(levels(p$data$SampleID), colnames(tse)[col_order])
    # Only rows sorted
    p <- plotHeatmap(tse, assay.type = "counts", sort.rows = "PCA_angle")
    expect_identical(levels(p$data$FeatureID), rownames(tse)[row_order])
    expect_identical(levels(p$data$SampleID), colnames(tse))
    # Disable column sorting
    p <- plotHeatmap(
        tse, assay.type = "counts", sort.method = "PCA_angle",
        sort.cols = NULL)
    expect_identical(levels(p$data$FeatureID), rownames(tse)[row_order])
    expect_identical(levels(p$data$SampleID), colnames(tse))
    # Tree defines row order
    expect_s3_class(
        plotHeatmap(
            tse, assay.type = "counts", show.tree = TRUE,
            sort.cols = "Heatmap"),
        "patchwork")
    expect_error(
        plotHeatmap(
            tse, assay.type = "counts", show.tree = TRUE,
            sort.rows = "PCA_angle"),
        "show.tree")
})
