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
    p <- plotForest(df, err.var = "err", df = 10, label.by = "CI")
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
        plotForest(df, err.var = "err", df = 10, ci.lower.var = "something"),
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
test_that("plotForest degrees of freedom for CIs from standard errors", {
    df <- data.frame(
        effect = c(1, -0.5, 0.2),
        err = c(0.2, 0.3, 0.1),
        dof = c(10, 25, 0),
        pval = c(0.01, 0.2, 0.5),
        row.names = c("a", "b", "c")
    )
    ci <- function(p) {
        d <- p[[1]]$data
        d[order(d$rownames), c("lower", "upper")]
    }
    # Per-row degrees of freedom; non-positive df gives no interval
    p <- plotForest(df, err.var = "err", df.var = "dof")
    res <- ci(p)
    q <- qt(0.975, c(10, 25))
    expect_equal(res$lower[1:2], df$effect[1:2] - q * df$err[1:2])
    expect_equal(res$upper[1:2], df$effect[1:2] + q * df$err[1:2])
    expect_true(all(is.na(res[3, ])))
    # A single value for all rows
    p <- plotForest(df, err.var = "err", df = 20)
    expect_equal(ci(p)$upper, df$effect + qt(0.975, 20) * df$err)
    # Without degrees of freedom: normal distribution with a warning
    expect_warning(
        p <- plotForest(df, err.var = "err"), "normal distribution")
    expect_equal(ci(p)$upper, df$effect + qnorm(0.975) * df$err)
    # Invalid input
    expect_error(plotForest(df, err.var = "err", df.var = "dof", df = 5))
    expect_error(plotForest(df, err.var = "err", df.var = "missing"))
    expect_error(plotForest(df, err.var = "err", df = -1))
    # SE method passes the degrees of freedom on
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, 3, 2,
            dimnames = list(rownames(df), c("s1", "s2")))),
        rowData = df)
    expect_s3_class(plotForest(se, err.var = "err", df = 20), "ggplot")
})

test_that("plotForest labels use fixed decimals and show small p-values", {
    df <- data.frame(
        effect = c(1.5, -0.45), lower = c(1.2, -1.594), upper = c(1.8, 0.7),
        pval = c(1e-6, 0.4281), row.names = c("a", "b"))
    p <- plotForest(df, label.by = c("CI", "P-Value"))
    # Collect the text of all label panels
    labs <- unlist(lapply(seq_along(p), function(i){
        q <- p[[i]]
        if( inherits(q, "patchwork") ){
            unlist(lapply(seq_along(q), function(j)
                tryCatch(layer_data(q[[j]])$label, error = function(e) NULL)))
        } else{
            tryCatch(layer_data(q)$label, error = function(e) NULL)
        }
    }))
    expect_true(all(c("1.50 (1.20, 1.80)", "-0.45 (-1.59, 0.70)") %in% labs))
    expect_true(all(c("<0.001*", "0.428") %in% labs))
})
