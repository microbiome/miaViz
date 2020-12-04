
context("plot tree")
test_that("plot tree", {
    # .check_tree_plot_switches
    expect_error(miaViz:::.check_tree_plot_switches(),
                 'argument "relabel_tree" is missing')
    expect_error(miaViz:::.check_tree_plot_switches(TRUE),
                 'argument "show_label" is missing')
    expect_error(miaViz:::.check_tree_plot_switches(TRUE, TRUE),
                 'argument "add_legend" is missing')
    expect_null(miaViz:::.check_tree_plot_switches(TRUE, TRUE, TRUE))
    expect_error(miaViz:::.check_tree_plot_switches("TRUE", TRUE, TRUE),
                 "'relabel_tree' must be either TRUE or FALSE")
    expect_error(miaViz:::.check_tree_plot_switches(TRUE, "TRUE", TRUE),
                 "'show_label' must be either TRUE or FALSE or named logical vector")
    expect_error(miaViz:::.check_tree_plot_switches(TRUE, c(TRUE,FALSE), TRUE),
                 "'show_label' must be either TRUE or FALSE or named logical vector")
    expect_error(miaViz:::.check_tree_plot_switches(TRUE, TRUE, "TRUE"),
                 "'add_legend' must be either TRUE or FALSE")
    #
    data("GlobalPatterns")
    x <- GlobalPatterns
    # .get_trimed_object_and_tree
    expect_error(miaViz:::.get_trimed_object_and_tree(),
                 'argument "object" is missing')
    actual <- miaViz:::.get_trimed_object_and_tree(x["549322",])
    expect_s3_class(actual$tree,"phylo")
    expect_s4_class(actual$object,"TreeSummarizedExperiment")
    expect_equal(unique(actual$tree$tip.label), c("549322"))
    actual <- miaViz:::.get_trimed_object_and_tree(x)
    expect_equal(actual$tree$tip.label, rownames(x))
    actual <- miaViz:::.get_trimed_object_and_tree(x, relabel = TRUE)
    expect_equal(actual$tree$tip.label[1L], "Class:Thermoprotei")
    #
    library(scater)
    library(mia)
    data(GlobalPatterns)
    altExp(GlobalPatterns,"genus") <- agglomerateByRank(GlobalPatterns,"Genus")
    altExp(GlobalPatterns,"genus") <- addPerFeatureQC(altExp(GlobalPatterns,"genus"))
    rowData(altExp(GlobalPatterns,"genus"))$log_mean <- log(rowData(altExp(GlobalPatterns,"genus"))$mean)
    top_taxa <- getTopTaxa(altExp(GlobalPatterns,"genus"),
                           method="mean",
                           top=100L,
                           abund_values="counts")
    #
    plot <- plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                tip_colour_by = "log_mean",
                tip_size_by = "detected")
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    # plot with tip labels
    plot <- plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                tip_colour_by = "log_mean",
                show_label = TRUE)
    expect_true(all(c("colour_by") %in% colnames(plot$data)))
    # plot with selected labels
    labels <- c("Genus:Providencia" = TRUE, "Genus:Morganella" = FALSE,
                "0.961.60" = TRUE)
    plot <- plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                tip_colour_by = "log_mean",
                tip_size_by = "detected",
                show_label = labels,
                layout="rectangular")
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    expect_equal(unique(plot$data$label), c("","Genus:Providencia","0.961.60"))
})
