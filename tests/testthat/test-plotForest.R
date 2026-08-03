test_that("plotForest", {
    # Make toy data
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    # Make toy stats
    est <- rnorm(nrow(tse), mean = 0)
    conf.lower <- est - 0.1 * rpois(length(est), lambda = 3)
    conf.upper <- est + 0.1 * rpois(length(est), lambda = 3)
    err <- 0.1 * rpois(length(est), lambda = 3)
    p.val <- runif(est)
    extra <- sample(100, length(est))
    # Compile into df
    df <- data.frame(
        effect = est,
        lower = conf.lower,
        upper = conf.upper,
        err = err,
        pval = p.val,
        extra = extra
    )
    # Add stats to side information
    rowData(tse) <- cbind(rowData(tse), df)
    colData(tse) <- cbind(colData(tse), df)
    # Check feature-wise behaviour (tree + forest + 3 text cols)
    p <- plotForest(tse, label.by = c("CI", "P-Value", "extra"))
    expect_s3_class(p, "ggplot")
    # Check sample-wise behaviour
    p <- plotForest(tse, by = 2, colour.by = "group")
    expect_s3_class(p, "ggplot")
    # Check df method
    p <- plotForest(df, colour.by = "extra")
    expect_s3_class(p, "ggplot")
    # Check TreeSE method with many args
    p <- plotForest(
        tse,
        label.by = c("CI", "P-Value", "extra"),
        colour.by = "var1",
        edge.colour.by = "var2"
    )
    expect_s3_class(p, "ggplot")
    # Reduce TreeSE to SE
    se <- SummarizedExperiment(
        assays = assays(tse),
        rowData = rowData(tse),
        colData = colData(tse)
    )
    # Check SE method
    p <- plotForest(se)
    expect_s3_class(p, "ggplot")
    # Check no errorbar when CI not defined
    p <- plotForest(
        se,
        effect.var = "effect",
        ci.lower.var = NULL
    )
    expect_s3_class(p, "ggplot")
    # Check paired forest plot (2 lines per row)
    df <- as.data.frame(colData(tse))
    # Compute CI from standard error
    p <- plotForest(df, err.var = "err", label.by = "CI")
    expect_s3_class(p, "ggplot")
    # Combine df for grouped analysis
    df2 <- rbind(df, df)
    df2$Group <- c(rep("A", nrow(df2) / 2), rep("B", nrow(df2) / 2))
    # Without colour and facet
    p_raw <- plotForest(df2, id.var = "ID")
    expect_length(unique(ggplot_build(p_raw)$data[[3]]$colour), 1L)
    expect_length(unique(ggplot_build(p_raw)$data[[3]]$PANEL), 1L)
    # With colour
    p_colour <- plotForest(df2, id.var = "ID", colour.by = "Group")
    expect_length(unique(ggplot_build(p_colour)$data[[3]]$colour), 2L)
    # With facet
    p_facet <- plotForest(df2, id.var = "ID", facet.by = "Group")
    expect_length(unique(ggplot_build(p_facet)$data[[3]]$PANEL), 2L)
    # Check argument values
    expect_error(
        plotForest(tse, tree.name = "wrong"),
        "'tree.name' must specify a tree from 'rowTreeNames(x)'.",
        fixed = TRUE
    )
    expect_error(
        plotForest(tse, show.tree = "not_a_logical"),
        "'show.tree' must be TRUE or FALSE.",
        fixed = TRUE
    )
    colTree(tse) <- NULL
    expect_error(
        plotForest(tse, by = 2L, show.tree = TRUE),
        "'tree.name' must specify a tree from 'colTreeNames(x)'.",
        fixed = TRUE
    )
    expect_warning(
        plotForest(tse, order.by = "effect", show.tree = TRUE),
        "'show.tree' is ignored when 'order.by' is defined.",
        fixed = TRUE
    )
    expect_error(
        plotForest(tse, layout = "circular", branch.length = "none"),
        "'layout', 'levels.rm' and 'branch.length' cannot be modified for this plot.",
        fixed = TRUE
    )
    expect_warning(
        plotForest(df, ci.lower.var = "wrong"),
        "'ci.lower.var' and 'ci.upper.var' are ignored when not found in x.",
        fixed = TRUE
    )
    expect_error(
        plotForest(df, ci.lower.var = NULL, label.by = "CI"),
        "To show the 'CI' label, x must include the variables specified with 'effect.var', 'ci.lower.var' and 'ci.upper.var'.",
        fixed = TRUE
    )
    expect_error(
        plotForest(df, pval.var = NULL, label.by = "P-Value"),
        "To show the 'P-Value' label, x must include the variable specified with 'pval.var'.",
        fixed = TRUE
    )
    expect_warning(
        plotForest(df, err.var = "err", ci.lower.var = "something"),
        "'ci.lower.var' and 'ci.upper.var' are ignored when 'err.var' is defined.",
        fixed = TRUE
    )
    expect_warning(
        plotForest(df, err.var = "wrong"),
        "'err.var' is ignored when not found in x.",
        fixed = TRUE
    )
    expect_error(
        plotForest(df, id.var = "wrong"),
        "'id.var' must be 'rownames' or a variable in x.",
        fixed = TRUE
    )
})