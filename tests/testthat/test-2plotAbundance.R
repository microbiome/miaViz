
context("plot abundance")
test_that("plot abundance", {
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .check_tree_plot_switches
    expect_error(miaViz:::.get_abundance_data(group = "Phylum"),
                 'argument "x" is missing')
    expect_error(miaViz:::.get_abundance_data(x, group = "Phylum"),
                 'argument "assay.type" is missing')
    expect_error(miaViz:::.get_abundance_data(x))
    actual <- miaViz:::.get_abundance_data(x, group = "Phylum", assay.type = "counts")
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
        x, group = "Phylum", assay.type = "counts", order.row.by = "abund",
        as.relative = FALSE)
    expect_true(max(actual3$Y) > 1)
    # .order_abundance_cols
    expect_error(miaViz:::.order_abundance_cols("meep"))
    expect_error(miaViz:::.order_abundance_cols("meep","meep2",x))
    # .get_feature_data
    expect_true(is.null(miaViz:::.get_feature_data()))
    expect_true(is.null(miaViz:::.get_feature_data(x)))
    expect_true(is.null(miaViz:::.get_feature_data(x,"meep")))
    actual <- miaViz:::.get_feature_data(x,"Primer")
    expect_true(is.list(actual))
    expect_named(actual,c("name","value"))
    # .get_features_data
    expect_error(miaViz:::.get_features_data(),
                 'argument "order_sample_by" is missing')
    expect_error(miaViz:::.get_features_data("Primer"),
                 'argument "order_sample_by" is missing')
    actual <- miaViz:::.get_features_data("Primer","SampleType",x)
    expect_true(is.data.frame(actual))
    expect_named(actual,c("SampleType","Primer"))
    # .order_abund_feature_data
    abund_data <- miaViz:::.get_abundance_data(
        x, group = "Phylum", assay.type = "counts")
    features_data <- miaViz:::.get_features_data("Primer","SampleType",x)
    expect_error(miaViz:::.order_abund_feature_data(abund_data,features_data),
                 'argument "order_sample_by" is missing')
    actual <- miaViz:::.order_abund_feature_data(abund_data,features_data,
                                                 "Primer")
    expect_true(is.list(actual))
    expect_named(actual,c("abund_data","features_data"))
    expect_equal(abund_data,actual$abund_data)
    expect_equal(features_data,actual$features_data)
    actual <- miaViz:::.order_abund_feature_data(abund_data,features_data,
                                                 "Primer",decreasing = FALSE)
    expect_error(expect_equal(abund_data,actual$abund_data))
    expect_error(expect_equal(features_data,actual$features_data))
    #
    expect_error(plot <- plotAbundance(x, assay.type="counts"))
    plot <- plotAbundance(x[1:20, ], assay.type="counts")
    expect_s3_class(plot,"ggplot")
    expect_named(plot$data,c("colour_by","X","Y"))
    plot <- plotAbundance(x, assay.type="counts", group = "Phylum",
                          col.var = "SampleType",
                          order.col.by = "SampleType")
    expect_true(is.list(plot))
    expect_s3_class(plot[[1]],"ggplot")
    #
    rowData(x)$Salame <- sample(letters[1:5], nrow(x), replace=TRUE)
    plot <- plotAbundance(x, assay.type="counts", group = "Salame",
                          col.var = "SampleType",
                          order.col.by = "SampleType")
    expect_true(is.list(plot))
    expect_s3_class(plot[[1]],"ggplot")
})
