
context("plot tree")
test_that("plot tree", {
    #
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .get_object_and_trimmed_tree
    expect_error(miaViz:::.get_object_and_trimmed_tree())
    actual <- miaViz:::.get_object_and_trimmed_tree(x["549322",])
    expect_s3_class(rowTree(actual),"phylo")
    expect_s4_class(actual,"TreeSummarizedExperiment")
    expect_equal(unique(rowTree(actual)$tip.label), c("549322"))
    actual <- miaViz:::.get_object_and_trimmed_tree(x)
    expect_equal(rowTree(actual)$tip.label, rownames(x))
    actual <- miaViz:::.get_object_and_trimmed_tree(x, relabel = TRUE)
    expect_equal(rowTree(actual)$tip.label[1L], "Class:Thermoprotei")
    #
    library(scater)
    library(mia)
    data(GlobalPatterns)
    altExp(GlobalPatterns,"genus") <- agglomerateByRank(GlobalPatterns,"Genus", make_unique = FALSE)
    altExp(GlobalPatterns,"genus") <- addPerFeatureQC(altExp(GlobalPatterns,"genus"))
    rowData(altExp(GlobalPatterns,"genus"))$log_mean <- log(rowData(altExp(GlobalPatterns,"genus"))$mean)
    top_taxa <- getTop(altExp(GlobalPatterns,"genus"),
                           method="mean",
                           top=100L,
                           assay.type="counts")
    #
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected"))
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    # plot with tip labels
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       show.label = TRUE))
    expect_true(all(c("colour_by") %in% colnames(plot$data)))
    # plot with selected labels
    labels <- c("Genus:Providencia", "Genus:Morganella", "0.961.60")
    plot <- expect_warning(plotRowTree(altExp(GlobalPatterns,"genus")[top_taxa,],
                                       tip.colour.by = "log_mean",
                                       tip.size.by = "detected",
                                       show.label = labels,
                                       layout="rectangular"))
    expect_true(all(c("colour_by", "size_by") %in% colnames(plot$data)))
    # Test that error occurs if tree.name is wrong
    expect_error( plotRowTree(GlobalPatterns, tree.name = "test") )
    expect_error( plotRowTree(GlobalPatterns, tree.name = NULL) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = c("test", "phylo")) )
    expect_error( plotColTree(GlobalPatterns, tree.name = 1) )
    expect_error( plotRowTree(GlobalPatterns, tree.name = TRUE) )
})
