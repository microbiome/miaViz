
context("plot graph")
test_that("plot graph", {
    library(tidygraph)
    # .check_graph_plot_switches
    expect_error(miaViz:::.check_graph_plot_switches(),
                 'argument "show.label" is missing')
    expect_error(miaViz:::.check_graph_plot_switches(TRUE),
                 'argument "add.legend" is missing')
    expect_error(miaViz:::.check_graph_plot_switches(TRUE, 1),
                 "'add_legend' must be either TRUE or FALSE")
    expect_null(miaViz:::.check_graph_plot_switches(1, FALSE))
    # .norm_layout_edge_type
    expect_error(miaViz:::.norm_layout_edge_type(),
                 'argument "edge.type" is missing')
    expect_error(miaViz:::.norm_layout_edge_type("meep","meep"),
                 "'arg' should be one")
    #
    data(GlobalPatterns)
    data("row_graph")
    data("col_graph")
    se <- GlobalPatterns
    # .get_graph_data
    actual <- miaViz:::.get_graph_data(row_graph)
    expect_s3_class(actual,c("tbl_graph","igraph"))
    actual2 <- miaViz:::.add_graph_node_labels(actual,TRUE)
    expect_s3_class(actual2$df,c("tbl_graph","igraph"))
    df <- actual2$df %>% activate("nodes") %>% as.data.frame()
    expect_named(df, c("name","label"))
    expect_equal(df$name,df$label)
    actual2 <- miaViz:::.add_graph_node_labels(actual,FALSE)
    expect_s3_class(actual2$df,c("tbl_graph","igraph"))
    df <- actual2$df %>% activate("nodes") %>% as.data.frame()
    expect_named(df, c("name","label"))
    expect_equal(df$name,df$label)
    show_label <- c(TRUE,
                    rep(FALSE,
                        nrow(actual %>% activate("nodes") %>% as.data.frame())))
    expect_error(miaViz:::.add_graph_node_labels(actual,show_label),
                 "If 'show.label' is logical")
    show_label <- show_label[-length(show_label)]
    actual2 <- miaViz:::.add_graph_node_labels(actual,show_label)
    expect_named(actual2$df %>% activate("nodes") %>% as.data.frame(),
                 c("name","label"))
    df <- actual2$df %>% activate("nodes") %>% as.data.frame()
    expect_named(df, c("name","label"))
    expect_true(all(is.na(df$label[-1])))
    # .colnames_tbl_graph
    expect_equal(miaViz:::.colnames_tbl_graph(actual,"nodes"),"name")
    expect_equal(miaViz:::.colnames_tbl_graph(actual,"edges"),
                 c("from","to","weight"))
    #
    genus <- agglomerateByRank(GlobalPatterns,"Genus",na.rm=TRUE)
    plot <- plotColGraph(col_graph,
                         genus,
                         colour.by = "SampleType",
                         edge.colour.by = "weight",
                         edge.width.by = "weight",
                         show.label = TRUE)
    expect_s3_class(plot,"ggplot")
    metadata(genus)$col_graph <- col_graph
    plot2 <- plotColGraph(genus,
                          name = "col_graph",
                         colour.by = "SampleType",
                         edge.colour.by = "weight",
                         edge.width.by = "weight",
                         show.label = TRUE)
    expect_true(all(plot$data == plot2$data))
})

test_that("colour.by works together with shape.by", {
    # Shapes from shape.by cannot be filled, so colour must be used
    args <- miaViz:::.get_point_args("colour_by", "shape_by", NULL)
    expect_false(args$fill)
    expect_true("colour" %in% names(args$args$mapping))
    expect_false("fill" %in% names(args$args$mapping))
    # Without shape.by, the default shape 21 is filled
    args <- miaViz:::.get_point_args("colour_by", NULL, NULL)
    expect_true(args$fill)
    expect_true("fill" %in% names(args$args$mapping))
    # Graph with 6 features in 2 groups
    tse <- TreeSummarizedExperiment(
        assays = list(counts = matrix(1, 6, 2,
            dimnames = list(paste0("f", 1:6), c("s1", "s2")))),
        rowData = DataFrame(group = rep(c("A", "B"), each = 3),
            type = rep(c("x", "y"), 3)))
    graph <- tidygraph::tbl_graph(
        nodes = data.frame(name = rownames(tse)),
        edges = data.frame(from = 1:5, to = 2:6), directed = FALSE)
    p <- plotRowGraph(graph, tse, colour.by = "group", shape.by = "type")
    node_layer <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)
    expect_true(all(c("colour", "shape") %in%
        names(node_layer[[1]]$mapping)))
    # Colours are not all the same
    expect_gt(length(unique(ggplot_build(p)$data[[2]]$colour)), 1)
})
