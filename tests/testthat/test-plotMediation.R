test_that("plot mediation", {
    # Load data from miaTime package
    data(hitchip1006, package = "miaTime")
    tse <- hitchip1006
    # Prepare data
    tse$bmi_group <- as.numeric(tse$bmi_group)
    tse <- agglomerateByRank(tse, rank = "Phylum")
    # Run mediation analysis
    tse <- addMediation(tse, outcome = "bmi_group", treatment = "nationality",
        mediator = "diversity", covariates = c("sex", "age"),
        treat.value = "Scandinavia", control.value = "CentralEurope",
        boot = TRUE, sims = 100)

    ### 1). TEST error messages ###
    expect_error(plotMediation(tse, "wrong_name"),
        "'name' must be in names(metadata(x)).", fixed = TRUE)
    expect_error(plotMediation(tse, layout = "barplot"),
        "'layout' must be either 'heatmap' or 'forest'.", fixed = TRUE)
    expect_error(plotMediation(tse, add.significance = "0.01"))

    ### 2). TEST plot layers ###
    df <- metadata(tse)[["mediation"]]
    p_heatmap <- plotMediation(df, layout = "heatmap")
    p_forest <- plotMediation(df, layout = "forest")
    heatmap_data <- ggplot_build(p_heatmap)[["data"]]
    forest_data <- ggplot_build(p_forest)[["data"]]
    # Check heatmap
    expect_equal(length(heatmap_data), 2)
    expect_equal(length(forest_data), 3)
    expect_true(all(grepl("^-?[0-9]+\\.[0-9]{2}\\*\\*\\*$",
        heatmap_data[[2]][["label"]])))
    # Without estimates only the symbols are shown
    p <- plotMediation(df, add.estimate = FALSE)
    expect_true(all(layer_data(p, 2L)[["label"]] == "***"))
    expect_error(plotMediation(df, add.estimate = "yes"),
        "'add.estimate' must be TRUE or FALSE.", fixed = TRUE)
    # Check forest plot
    expect_true(all(forest_data[[2]][["xmin"]] ==
        df[c("acme_lower", "ade_lower", "total_lower")]))
    expect_true(all(forest_data[[2]][["xmax"]] ==
        df[c("acme_upper", "ade_upper", "total_upper")]))
    expect_true(all(forest_data[[2]][["x"]] ==
        df[c("acme", "ade", "total")]))
})

test_that("plotMediation forest shows components in their natural order", {
    # Mediation table in the format of mia, sorted by significance
    .med_df <- function(mediator){
        n <- length(mediator)
        df <- data.frame(mediator = mediator)
        for( col in c("acme", "ade", "total") ){
            df[[col]] <- seq_len(n)
            df[[paste0(col, "_lower")]] <- seq_len(n) - 1
            df[[paste0(col, "_upper")]] <- seq_len(n) + 1
            df[[paste0(col, "_pval")]] <- seq_len(n) / 100
            df[[paste0(col, "_padj")]] <- seq_len(n) / 100
        }
        df
    }
    .facets <- function(p){
        as.character(ggplot_build(p)$layout$layout[["mediator"]])
    }
    p <- plotMediation(.med_df(c("PCA10", "PCA1", "PCA3", "PCA2")),
        layout = "forest")
    expect_equal(.facets(p), c("PCA1", "PCA2", "PCA3", "PCA10"))
    # Other mediators keep the order of significance
    p <- plotMediation(.med_df(c("c", "a", "b")), layout = "forest")
    expect_equal(.facets(p), c("c", "a", "b"))
})

test_that("plotMediation heatmap shares one symmetric scale", {
    df <- data.frame(mediator = c("a", "b"))
    vals <- list(acme = c(-0.004, 0.01), ade = c(-0.2, -0.1),
        total = c(-0.3, -0.3))
    for( col in names(vals) ){
        df[[col]] <- vals[[col]]
        df[[paste0(col, "_lower")]] <- vals[[col]] - 0.1
        df[[paste0(col, "_upper")]] <- vals[[col]] + 0.1
        df[[paste0(col, "_pval")]] <- c(0.01, 0.5)
        df[[paste0(col, "_padj")]] <- c(0.01, 0.5)
    }
    p <- plotMediation(df)
    expect_equal(p$scales$get_scales("fill")$get_limits(), c(-0.3, 0.3))
    labs <- layer_data(p, 2L)[["label"]]
    # No negative zero, and the symbol follows the estimate
    expect_true("0.00*" %in% labs)
    expect_false(any(grepl("^-0\\.00", labs)))
    expect_equal(p$labels$fill, "Effect")
    # Numeric p-values go on their own line
    p <- plotMediation(df, add.significance = "numeric")
    expect_true("-0.30\np=0.010" %in% layer_data(p, 2L)[["label"]])
})
