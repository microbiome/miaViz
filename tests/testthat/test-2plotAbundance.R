
context("plot abundance")
test_that("plot abundance", {
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .check_tree_plot_switches
    expect_error(
        miaViz:::.get_abundance_data(group = "Phylum"),
        'argument "x" is missing')
    expect_error(
        miaViz:::.get_abundance_data(x, group = "Phylum"),
        'argument "assay.type" is missing')
    expect_error(miaViz:::.get_abundance_data(x))
    actual <- miaViz:::.get_abundance_data(
        x, group = "Phylum", assay.type = "counts")
    expect_s3_class(actual,"tbl_df")
    expect_true( all(c("colour_by","X","Y") %in% colnames(actual)) )
    expect_true(is.factor(actual$colour_by))
    expect_true(is.factor(actual$X))
    expect_true(is.numeric(actual$Y))
    expect_error(miaViz:::.order_abundance_rows(actual, order.row.by = "meep"))
    actual2 <- miaViz:::.order_abundance_rows(actual, order.row.by = "abund")
    expect_equal(as.character(actual[1,1,drop=TRUE]),"ABY1_OD1")
    expect_equal(levels(actual2[["colour_by"]])[[1]],"Proteobacteria")
    actual3 <- miaViz:::.get_abundance_data(
        x, group = "Phylum", assay.type = "counts",
        order.row.by = "abund",
        as.relative = FALSE)
    expect_true(max(actual3$Y, na.rm = TRUE) > 1)
    # .order_abundance_cols
    expect_error(miaViz:::.order_abundance_cols("meep"))
    expect_error(miaViz:::.order_abundance_cols("meep","meep2",x))
    #
    expect_error(plot <- plotAbundance(x, assay.type="counts"))
    plot <- plotAbundance(x[1:20, ], assay.type="counts")
    expect_s3_class(plot,"ggplot")
    expect_true(all(c("colour_by","X","Y") %in% colnames(plot$data)))
    plot <- plotAbundance(
        x, assay.type = "counts", group = "Phylum", col.var = "SampleType",
        order.col.by = "SampleType")
    expect_true(is.list(plot))
    expect_s3_class(plot[[1]],"ggplot")
    # Check that the grouping is correct
    expect_true(all(levels(plot$data$colour_by) %in% unique(rowData(x)$Phylum)))
    #
    rowData(x)$Salame <- sample(letters[1:5], nrow(x), replace=TRUE)
    plot <- plotAbundance(
        x, assay.type="counts", group = "Salame", col.var = "SampleType",
        order.col.by = "SampleType")
    expect_true(is.list(plot))
    expect_s3_class(plot[[1]],"ggplot")
    # Expect error since too many rows
    expect_error(plotAbundance(x, assay.type="counts"))
    # Mock data for paired sample testing
    df <- data.frame(
        sample = c("S1", "S2", "S3"),
        time_point = c(1, 1, 2),
        value = c(10, 15, 20)
    )
    # Test paired sample logic with correct input
    paired_df <- miaViz:::.add_paired_samples(df, paired = TRUE, order.col.by = "sample")
    expect_equal(paired_df, df)  # Should still have 3 rows (no missing)
    # Test with invalid `paired` input
    expect_error(miaViz:::.add_paired_samples(df, paired = "INVALID"))
    # Test edge case with empty data
    empty_df <- data.frame()
    expect_equal(nrow(miaViz:::.add_paired_samples(empty_df)), 0)
    #
    # Test with correct input
    actual <- miaViz:::.get_abundance_data(
        x, group = "Phylum", assay.type = "counts")
    expect_s3_class(actual, "tbl_df")
    expect_true(all(c("colour_by", "X", "Y") %in% colnames(actual)))
    # Test error when missing necessary arguments
    expect_error(miaViz:::.get_abundance_data(x[1:10, ], group = "Phylum"),
        "argument \"assay.type\" is missing, with no default")
    expect_error(miaViz:::.get_abundance_data(x[1:10, ]),
        "argument \"assay.type\" is missing, with no default")
    # Edge case with empty dataset
    empty_x <- x[NULL, ]  # Empty dataset
    expect_error(miaViz:::.get_abundance_data(
        empty_x, group = "Phylum", assay.type = "counts"))
    # Test with correct input
    df <- .get_abundance_data(x[1:20, ], "counts")
    actual <- miaViz:::.order_abundance_cols(
        df, order.col.by = "SampleType")
    expect_s3_class(actual, "tbl_df")
    # Test error for incorrect ordering variable
    expect_error(miaViz:::.order_abundance_cols(
        df, order.col.by = "InvalidColumn"))
    # Test edge case with empty dataset
    empty_df <- df[NULL, ]  # Empty dataset
    expect_equal(nrow(miaViz:::.order_abundance_cols(
        empty_df, order.col.by = "SampleType")), 0)
    # Edge case: Empty dataset
    x_empty <- x[NULL, ]
    expect_error(plotAbundance(x_empty, assay.type = "counts"))
    # Edge case: Dataset missing the assay.type column
    x_missing <- x
    assays(x_missing)$counts <- NULL  # Remove assay.type column
    expect_error(plotAbundance(x_missing[1:10, ], assay.type = "counts"))
    # Test with numeric data as grouping factor
    x_numeric_group <- x[1:20, ]
    rowData(x_numeric_group)$Group <- 1:20  # Numeric group instead of factor
    expect_error(plot <- plotAbundance(
        x_numeric_group, assay.type="counts", group = "Group"))
})
