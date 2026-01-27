test_that("plotForest", {
    # Make toy data
    tse <- makeTSE(nrow = 10L, ncol = 10L)
    # Make toy stats
    est <- rnorm(nrow(tse), mean = 0)
    conf.lower <- est - 0.1 * rpois(length(est), lambda = 3)
    conf.upper <- est + 0.1 * rpois(length(est), lambda = 3)
    p.val <- runif(est)
    extra <- sample(100, length(est))
    # Compile into df
    df <- data.frame(
        effect = est,
        lower = conf.lower,
        upper = conf.upper,
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
        show.tree = FALSE,
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
    df2 <- rbind(df, df)
    df2$Group <- c(rep("A", nrow(df2) / 2), rep("B", nrow(df2) / 2))
    # With colour
    p <- plotForest(df2, id.var = "ID", colour.by = "Group")
    expect_s3_class(p, "ggplot")
    # Without colour
    p <- plotForest(df2, id.var = "ID")
    expect_s3_class(p, "ggplot")
})