
context("plot graph")
test_that("plot graph", {
    library(tidygraph)
    # .check_graph_plot_switches
    expect_error(miaViz:::.check_graph_plot_switches(),
                 'argument "show.label" is missing')
    expect_error(miaViz:::.check_graph_plot_switches(TRUE),
                 'argument "add.legend" is missing')
    expect_error(miaViz:::.check_graph_plot_switches(TRUE, 1),
                 "'add.legend' must be either TRUE or FALSE")
    expect_null(miaViz:::.check_graph_plot_switches(1, FALSE))
    # .norm_layout_edge_type
    expect_error(miaViz:::.norm_layout_edge_type(),
                 'argument "edge_type" is missing')
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
    show.label <- c(TRUE,
                    rep(FALSE,
                        nrow(actual %>% activate("nodes") %>% as.data.frame())))
    expect_error(miaViz:::.add_graph_node_labels(actual,show.label),
                 "If 'show.label' is logical")
    show.label <- show.label[-length(show.label)]
    actual2 <- miaViz:::.add_graph_node_labels(actual,show.label)
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
