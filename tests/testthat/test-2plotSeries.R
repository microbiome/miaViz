context("plot series")

test_that("plot series", {
    # Load data from miaTime package
    skip_if_not_installed("miaTime")
    data(SilvermanAGutData, package = "miaTime")
    tse <- SilvermanAGutData
    tse_sub <- tse[1:5, ]

    # Expect error
    expect_error(plotSeries())
    expect_error(plotSeries(tse_sub))

    # Expect output to be a ggplot object
    expect_s3_class(plotSeries(tse_sub, assay.type = "counts", time.col = "DAY_ORDER"), "ggplot")

    # Expect warning when over 10 taxa, expect error when over 20 taxa
    tse_sub <- tse[1:11, ]
    expect_warning(plotSeries(tse_sub, assay.type = "counts", time.col = "DAY_ORDER"),
                   "Over 10 taxa selected.")
    tse_sub <- tse[1:21, ]
    expect_error(plotSeries(tse_sub, assay.type = "counts", time.col = "DAY_ORDER"),
                 "Over 20 taxa selected. 20 or under allowed.")

    ##################### Test .get_series_data #################################
    tse_sub <- tse[1:10, ]
    expect_error(miaViz:::.get_series_data(tse_sub, "counts", NULL, "DAY_ORDER"),
                 "argument \"group.by\" is missing, with no default")

    ################### Test .series_plotter ###########################
    # Get data from colData
    x_test <- colData(tse)[, "DAY_ORDER"]

    # Get data from rowData
    colour_by_test <- as.character(rowData(tse)[, "Phylum"])
    linetype_by_test <- as.character(rowData(tse)[, "Family"])
    size_by_test <- as.character(rowData(tse)[, "Kingdom"])
    colour_linetype_and_size <- data.frame(colour_by = colour_by_test,
                                           linetype_by = linetype_by_test,
                                           size_by = size_by_test, stringsAsFactors = FALSE)

    # Get data from function
    data_from_function <- miaViz:::.get_series_data(
        tse, assay.type = "counts", col.var = NULL, time.col = "DAY_ORDER",
        colour.by = "Phylum", linetype.by = "Family", size.by = "Kingdom",
        group.by = NULL)

    # Extracting the relevant columns from the data
    series_data_test <- data_from_function$plot_data
    feature_data_test_without_rownames <- colour_linetype_and_size
    rownames(feature_data_test_without_rownames) <- NULL

    # Expect that data frames are equal in terms of structure
    expect_true("X" %in% names(series_data_test))
    expect_true("Y" %in% names(series_data_test))

    ################## Test .melt_series_data ##################################
    melted <- miaViz:::.get_series_data(
        tse_sub, "counts", col.var = NULL, group.by = NULL, "DAY_ORDER",
        "Phylum", "Family", "Kingdom")

    # Convert melted data to data frame
    melted <- as.data.frame(melted$plot_data)

    # Check 10 different random combinations
    for (i in 1:10) {
        # Get random taxa name
        taxa <- rownames(tse)[sample(length(rownames(tse)), 1)]
        # Get random sample name
        sample <- colnames(tse)[sample(length(colnames(tse)), 1)]
        # Get time point of the sample
        timepoint <- colData(tse)[sample, "DAY_ORDER"]

        # Ensure timepoint is of length 1
        if (length(timepoint) != 1) next

        # Get sample names from that time point
        sample_names <- rownames(colData(tse)[colData(tse)[, "DAY_ORDER"] == timepoint, ])

        # If there are no sample names for this timepoint, skip this iteration
        if (length(sample_names) == 0) next

        # Get values from assay data, and take the average
        assay_mean_value <- mean(assays(tse)$counts[taxa, sample_names], na.rm = TRUE)
        assay_sd_value <- sd(assays(tse)$counts[taxa, sample_names], na.rm = TRUE)

        # Get value from melted data
        melted_assay_sd_value <- as.double(melted[melted[, "feature"] == taxa & melted[, "X"] == timepoint, "sd"])
        melted_assay_mean_value <- as.double(melted[melted[, "feature"] == taxa & melted[, "X"] == timepoint, "Y"])

        # Ensure values are not NA and have length 1 before comparison
        if (all(!is.na(c(melted_assay_mean_value, assay_mean_value, melted_assay_sd_value, assay_sd_value))) &&
            length(melted_assay_mean_value) == 1 && length(assay_mean_value) == 1 &&
            length(melted_assay_sd_value) == 1 && length(assay_sd_value) == 1) {
            expect_equal(round(melted_assay_mean_value, 2), round(assay_mean_value, 2))
            expect_equal(round(melted_assay_sd_value, 2), round(assay_sd_value, 2))
        }

        # Phylum was chosen to be value of colour_by
        rowData_colour_by_value <- as.character(rowData(tse)[taxa, "Phylum"])
        melted_colour_by_value <- as.character(melted[melted[, "feature"] == taxa & melted[, "X"] == timepoint, "colour_by"])

        ######### Expect that colour_by values are melted correctly ########
        if (length(rowData_colour_by_value) == 1 && length(melted_colour_by_value) == 1 &&
            !is.na(rowData_colour_by_value) && !is.na(melted_colour_by_value)) {
            expect_equal(melted_colour_by_value, rowData_colour_by_value)
        }

        # Family was chosen to be value of linetype_by
        rowData_linetype_by_value <- as.character(rowData(tse)[taxa, "Family"])
        melted_linetype_by_value <- as.character(melted[melted[, "feature"] == taxa & melted[, "X"] == timepoint, "linetype_by"])

        ######### Expect that linetype_by is melted correctly ########
        if (length(rowData_linetype_by_value) == 1 && length(melted_linetype_by_value) == 1 &&
            !is.na(rowData_linetype_by_value) && !is.na(melted_linetype_by_value)) {
            expect_equal(melted_linetype_by_value, rowData_linetype_by_value)
        }

        # Kingdom was chosen to be value of size_by
        rowData_size_by_value <- as.character(rowData(tse)[taxa, "Kingdom"])
        melted_size_by_value <- as.character(melted[melted[, "feature"] == taxa & melted[, "X"] == timepoint, "size_by"])

        ######### Expect that size_by is melted correctly #########
        if (length(rowData_size_by_value) == 1 && length(melted_size_by_value) == 1 &&
            !is.na(rowData_size_by_value) && !is.na(melted_size_by_value)) {
            expect_equal(melted_size_by_value, rowData_size_by_value)
        }
    }

    # Check that the function works with col.var
    tse <- addAlpha(tse, index = "shannon")
    plotSeries(
        tse, assay.type = "counts", col.var = "shannon",
        time.col = "DAY_ORDER") |>
        expect_error()
    plotSeries(
        tse, colour.by = "Family", col.var = "shannon",
        time.col = "DAY_ORDER") |>
        expect_error()
    p <- plotSeries(tse, col.var = "shannon", time.col = "DAY_ORDER")
    expect_s3_class(p, "ggplot")
})
