
context("plot tree")
test_that("plot tree", {
    # .check_tree_plot_switches
    expect_error(miaViz:::.check_tree_plot_switches(),
                 'argument "layout" is missing')
    expect_error(miaViz:::.check_tree_plot_switches(TRUE),
                 "'layout' must be a single character value")
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE),
                 'argument "remove.levels" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE),
                 'argument "order.tree" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE),
                 'argument "show.label" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE),
                 'argument "show.highlights" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE,
                                                    TRUE),
                 'argument "show.highlight.label" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE,
                                                    TRUE, TRUE),
                 'argument "abbr.label" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE,
                                                    TRUE, TRUE,TRUE),
                 'argument "add.legend" is missing')
    expect_error(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE,
                                                    TRUE, TRUE,TRUE),
                 'argument "add.legend" is missing')
    expect_null(miaViz:::.check_tree_plot_switches("a",TRUE,TRUE,TRUE,TRUE,
                                                   TRUE, TRUE,TRUE, TRUE))
    #
    expect_error(miaViz:::.check_tree_plot_switches("A","TRUE", TRUE,TRUE,TRUE,
                                                    TRUE, TRUE,TRUE, TRUE),
                 "'relabel.tree' must be either TRUE or FALSE")
    expect_error(miaViz:::.check_tree_plot_switches("A",TRUE, 2, TRUE,TRUE,
                                                    TRUE, TRUE,TRUE, TRUE),
                 "'remove.levels' must be either TRUE or FALSE")
    expect_error(miaViz:::.check_tree_plot_switches("A",TRUE, TRUE, 2, TRUE,
                                                    TRUE, TRUE,TRUE, TRUE),
                 "'order.tree' must be either TRUE or FALSE")
    expect_null(miaViz:::.check_tree_plot_switches("A",TRUE, TRUE, TRUE, 2,
                                                    TRUE, TRUE,TRUE, TRUE))
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .get_object_and_trimmed_tree
    expect_error(miaViz:::.get_object_and_trimmed_tree(),
                 'argument "object" is missing')
    actual <- miaViz:::.get_object_and_trimmed_tree(x["549322",])
    expect_s3_class(actual$tree,"phylo")
    expect_s4_class(actual$object,"TreeSummarizedExperiment")
    expect_equal(unique(actual$tree$tip.label), c("549322"))
    actual <- miaViz:::.get_object_and_trimmed_tree(x)
    expect_equal(actual$tree$tip.label, rownames(x))
    actual <- miaViz:::.get_object_and_trimmed_tree(x, relabel = TRUE)
    expect_equal(actual$tree$tip.label[1L], "Class:Thermoprotei")
    #
    library(scater)
    library(mia)
    data(GlobalPatterns)
    altExp(GlobalPatterns,"genus") <- agglomerateByRank(GlobalPatterns,"Genus", make_unique = FALSE)
    altExp(GlobalPatterns,"genus") <- addPerFeatureQC(altExp(GlobalPatterns,"genus"))
    rowData(altExp(GlobalPatterns,"genus"))$log_mean <- log(rowData(altExp(GlobalPatterns,"genus"))$mean)
    top_taxa <- getTopFeatures(altExp(GlobalPatterns,"genus"),
                           method="mean",
                           top=100L,
                           assay.type="counts")
    #
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected"))
    expect_true(all(c("colour.by", "size.by") %in% colnames(plot$data)))
    # plot with tip labels
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       show.label = TRUE))
    expect_true(all(c("colour.by") %in% colnames(plot$data)))
    # plot with selected labels
    labels <- c("Genus:Providencia", "Genus:Morganella", "0.961.60")
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected",
                                       show.label = labels,
                                       layout="rectangular"))
    expect_true(all(c("colour.by", "size.by") %in% colnames(plot$data)))
    # Test that error occurs if tree.name is wrong
    expect_error( plotRowTree(GlobalPatterns, tree.name = "test") )
    expect_error( plotRowTree(GlobalPatterns, tree.name = NULL) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = c("test", "phylo")) )
    expect_error( plotColTree(GlobalPatterns, tree.name = 1) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = TRUE) )
})
