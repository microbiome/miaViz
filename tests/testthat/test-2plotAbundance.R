
context("plot abundance")
test_that("plot abundance", {
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .check_tree_plot_switches
    expect_error(miaViz:::.get_abundance_data(),
                 'argument "x" is missing')
    expect_error(miaViz:::.get_abundance_data(x),
                 'argument "assay.type" is missing')
    expect_error(miaViz:::.get_abundance_data(x,assay.type = "counts"),
                 'argument "rank" is missing')
    actual <- miaViz:::.get_abundance_data(x,"Phylum","counts")
    expect_s3_class(actual,"tbl_df")
    expect_named(actual,c("colour.by","X","Y"))
    expect_true(is.factor(actual$colour.by))
    expect_true(is.factor(actual$X))
    expect_true(is.numeric(actual$Y))
    expect_error(miaViz:::.get_abundance_data(x,"Phylum","counts",
                                              order.row.by = "meep"))
    actual2 <- miaViz:::.get_abundance_data(x,"Phylum","counts",
                                           order.row.by = "abund")
    expect_equal(as.character(actual[1,1,drop=TRUE]),"ABY1_OD1")
    expect_equal(as.character(actual2[1,1,drop=TRUE]),"Proteobacteria")
    actual3 <- miaViz:::.get_abundance_data(x,"Phylum","counts",
                                            order.row.by = "abund",
                                            use.relative = FALSE)
    expect_true(max(actual3$Y) > 1)
    # .norm_order_sample_by
    expect_true(is.null(miaViz:::.norm_order_sample_by(NULL)))
    expect_error(miaViz:::.norm_order_sample_by("meep"),
                 'argument "factors" is missing')
    expect_error(miaViz:::.norm_order_sample_by("meep","meep2",x),
                 "'order.col.by' must be a single non-empty character value")
    expect_equal(miaViz:::.norm_order_sample_by("Primer","meep2",x),"Primer")
    # .get_feature_data
    expect_true(is.null(miaViz:::.get_feature_data()))
    expect_true(is.null(miaViz:::.get_feature_data(x)))
    expect_true(is.null(miaViz:::.get_feature_data(x,"meep")))
    actual <- miaViz:::.get_feature_data(x,"Primer")
    expect_true(is.list(actual))
    expect_named(actual,c("name","value"))
    # .get_features_data
    expect_error(miaViz:::.get_features_data(),
                 'argument "order.col.by" is missing')
    expect_error(miaViz:::.get_features_data("Primer"),
                 'argument "order.col.by" is missing')
    actual <- miaViz:::.get_features_data("Primer","SampleType",x)
    expect_true(is.data.frame(actual))
    expect_named(actual,c("SampleType","Primer"))
    # .order_abund_feature_data
    abund_data <- miaViz:::.get_abundance_data(x,"Phylum","counts")
    features_data <- miaViz:::.get_features_data("Primer","SampleType",x)
    expect_error(miaViz:::.order_abund_feature_data(abund_data,features_data),
                 'argument "order.col.by" is missing')
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
    plot <- plotAbundance(x, assay.type="counts")
    expect_s3_class(plot,"ggplot")
    expect_named(plot$data,c("colour.by","X","Y"))
    plot <- plotAbundance(x, assay.type="counts", rank = "Phylum",
                          features = "SampleType",
                          order.col.by = "SampleType")
    expect_true(is.list(plot))
    expect_s3_class(plot[[1]],"ggplot")
})
