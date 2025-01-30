test_that("plot mediation", {
    # Load data from miaTime package
    skip_if_not_installed("miaTime")
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
    expect_true(all(heatmap_data[[2]][["label"]] == "***"))
    # Check forest plot
    expect_true(all(forest_data[[2]][["xmin"]] ==
        df[c("acme_lower", "ade_lower", "total_lower")]))
    expect_true(all(forest_data[[2]][["xmax"]] ==
        df[c("acme_upper", "ade_upper", "total_upper")]))
    expect_true(all(forest_data[[2]][["x"]] ==
        df[c("acme", "ade", "total")]))
})
