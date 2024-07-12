
context("tree data")
test_that("tree data", {
    data(GlobalPatterns)
    x <- GlobalPatterns
    # .get_tree_data
    expect_error(miaViz:::.get_tree_data(),
                 'argument "tree" is missing')
    actual <- miaViz:::.get_tree_data(rowTree(x))
    expect_s3_class(actual, c("data.frame","tbl_tree"))
    # .clean_tree_data
    actual$test <- "test"
    actual <- miaViz:::.clean_tree_data(actual)
    td <- actual
    expect_equal(colnames(actual), c("parent","node","branch.length","label"))
    # .norm_other_fields
    expect_error(miaViz:::.norm_other_fields(),
                 'argument "other.fields" is missing')
    expect_error(miaViz:::.norm_other_fields(list(a=c(1,2), b=c(2,3,4))))
    expect_null(miaViz:::.norm_other_fields(NULL))
    expect_null(miaViz:::.norm_other_fields(data.frame()))
    other_fields <- data.frame(a = "a", row.names = "549322")
    actual <- miaViz:::.norm_other_fields(other_fields)
    expect_equal(colnames(actual), c("node","a"))
    other_fields$node <- rownames(other_fields)
    rownames(other_fields) <- NULL
    expect_equal(actual, miaViz:::.norm_other_fields(other_fields))
    # .norm_id_col_of_other_fields
    expect_error(miaViz:::.norm_id_col_of_other_fields(),
                 'argument "other.fields" is missing')
    expect_error(miaViz:::.norm_id_col_of_other_fields(other_fields))
    expect_warning(miaViz:::.norm_id_col_of_other_fields(other_fields, td),
                   "Not all 'node' values found in tree data")
    other_fields$label <- other_fields$node
    other_fields$node <- NULL
    other_fields <- miaViz:::.norm_other_fields(other_fields)
    expect_equal(other_fields, miaViz:::.norm_id_col_of_other_fields(other_fields, td))
    expect_null(miaViz:::.norm_id_col_of_other_fields(NULL))
    # .combine_tree_data_and_other_fields
    expect_error(miaViz:::.combine_tree_data_and_other_fields(td),
                 'argument "other.fields" is missing')
    expect_equal(td, miaViz:::.combine_tree_data_and_other_fields(td, NULL))
    actual <- miaViz:::.combine_tree_data_and_other_fields(td, other_fields)
    expect_equal(colnames(actual),
                 c("parent","node","branch.length","label","a"))
    expect_equal(actual$a[1L],"a")
    expect_true(is.na(actual$a[2L]))
    # combineTreeData
    expect_equal(combineTreeData(rowTree(x)), tidytree::as.treedata(td))
    expect_equal(combineTreeData(rowTree(x), other_fields),
                 tidytree::as.treedata(actual))
    # Test situation when tree.name is wrong
    expect_equal( colTreeData(x, tree.name = "phylo"), NULL )
    expect_equal( rowTreeData(x, tree.name = "phylo.1"), NULL )
    expect_error( colTreeData(x, tree.name = 1) )
    expect_error( rowTreeData(x, tree.name = TRUE) )
    expect_error( rowTreeData(x, tree.name = c("test", "test2")) )
    td <- rowTreeData(x)
    rowTreeData(x, "test") <- td
    expect_equal(names(x@rowTree), c("phylo", "test"))
    data(esophagus)
    tse <- mia::mergeSEs(esophagus, GlobalPatterns)
    expect_equal( rowTreeData(tse, "phylo"), rowTreeData(esophagus, "phylo") )
})
