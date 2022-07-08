
context("plot series")
test_that("plot series", {
    
    tse <- microbiomeDataSets::SilvermanAGutData()
    tse_sub <- tse[1:5]
    
    # Expect error
    expect_error(plotSeries())
    expect_error(plotSeries(tse_sub))
    
    # Expect output
    expect_type(plotSeries(tse_sub, assay_name = "counts", x = "DAY_ORDER"), "list")
    
    # Expect warning when over 10 taxa, expect error when over 20 taxa
    tse_sub <- tse[1:11]
    expect_warning(plotSeries(tse_sub, assay_name = "counts", x = "DAY_ORDER"))
    tse_sub <- tse[1:21]
    expect_error(plotSeries(tse_sub, assay_name = "counts", x = "DAY_ORDER"))
    
    ##################### Test .get_assay_data #################################
    tse_sub <- tse[1:10]
    expect_equal(miaViz:::.get_assay_data(tse_sub, "counts"), assays(tse_sub)$counts)
    
    ################### Test .incorporate_series_vis ###########################
    # Get data from colData
    x_test <- colData(tse)[, "DAY_ORDER"]
    
    # Get data from rowData
    colour_by_test <- rowData(tse)[, "Phylum"]
    linetype_by_test <- rowData(tse)[,"Family"]
    size_by_test <- rowData(tse)[,"Kingdom"]
    colour_linetype_and_size <- data.frame(colour_by = colour_by_test, 
                                          linetype_by = linetype_by_test, 
                                          size_by = size_by_test)
    
    # Get data from function
    data_from_function <- miaViz:::.incorporate_series_vis(tse, x = "DAY_ORDER", 
                                                           colour_by = "Phylum", linetype_by = "Family",
                                                           size_by = "Kingdom")
    
    # Divide list to 2 data frames
    series_data_test <- data_from_function$series_data
    feature_data_test_without_rownames <- data_from_function$feature_data
    rownames(feature_data_test_without_rownames) <- NULL
    feature_data_test <- data_from_function$feature_data
    
    # Expect that data frames are equal
    expect_equal(series_data_test, x_test)
    expect_equal(feature_data_test_without_rownames, colour_linetype_and_size)
    
    ################## Test .melt_series_data ##################################
    melted <- as.data.frame(miaViz:::.melt_series_data(assays(tse)$counts, series_data_test, feature_data_test))
    
    # Check 10 different random combinations
    for ( i in c(1:10) ){
        # Get random taxa name
        taxa <- rownames(tse)[sample(length(rownames(tse)), 1)]
        # Get random sample name
        sample <- colnames(tse)[sample(length(colnames(tse)), 1)]
        # Get time point of the sample
        timepoint <- colData(tse)[rownames(colData(tse)) == sample, ]["DAY_ORDER"]
        
        # Get sample names from that time point
        sample_names <- rownames(colData(tse)[colData(tse)["DAY_ORDER"] == timepoint, ])
        
        # Get values from assay data, and take the average
        assay_mean_value = mean(assays(tse)$counts[taxa, sample_names], na.rm = TRUE)
        assay_sd_value = sd(assays(tse)$counts[taxa, sample_names], na.rm = TRUE)
        
        # Get value from melted data. Convert them to double, because they are characters
        melted_assay_sd_value <- as.double(melted[melted[, "feature"] == taxa, ][melted[melted[, "feature"] == taxa, ]["X"] == timepoint[,1]][3])
        melted_assay_mean_value <- as.double(melted[melted[, "feature"] == taxa, ][melted[melted[, "feature"] == taxa, ]["X"] == timepoint[,1]][4])
        
        ######### Expect that assay data is melted correctly ########
        expect_equal(round(melted_assay_mean_value, 2), round(assay_mean_value,2))
        expect_equal(round(melted_assay_sd_value, 2), round(assay_sd_value), 2)
        
        # Phylum was chosen to be value of colour_by
        rowData_colour_by_value <- rowData(tse)[taxa, "Phylum"]
        melted_colour_by_value <- melted[melted[, "feature"] == taxa, ][melted[melted[, "feature"] == taxa, ]["X"] == timepoint[,1]][5]
        
        ######### Expect that colour_by values are melted correctly ########
        expect_equal(melted_colour_by_value, rowData_colour_by_value)
        
        # Family was chosen to be value of linetype_by
        rowData_linetype_by_value <- rowData(tse)[taxa, "Family"]
        melted_linetype_by_value <- melted[melted[, "feature"] == taxa, ][melted[melted[, "feature"] == taxa, ]["X"] == timepoint[,1]][6]
        
        ######### Expect that linetype_by is melted correctly ########
        expect_equal(melted_linetype_by_value, rowData_linetype_by_value)
        
        # Kingdom was chosen to be value of size_by
        rowData_size_by_value <- rowData(tse)[taxa, "Kingdom"]
        melted_size_by_value <- melted[melted[, "feature"] == taxa, ][melted[melted[, "feature"] == taxa, ]["X"] == timepoint[,1]][7]
        
        ######### Expect that sizetype_by is melted correctly #########
        expect_equal(melted_size_by_value, rowData_size_by_value)
        
    }
    
})
