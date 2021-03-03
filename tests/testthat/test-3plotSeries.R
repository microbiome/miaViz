
context("plot series")
test_that("plot series", {
    
    data("GlobalPatterns")
    x <- GlobalPatterns
    x_sub <- x[1:5]
    
    # Expect error
    expect_error(plotSeries())
    expect_error(plotSeries(x_sub))
    expect_error(plotSeries(x_sub, abund_values = "counts"))
    expect_error(plotSeries(x_sub, X = "SampleType"))
    
    # Expect output
    expect_type(plotSeries(x_sub, abund_values = "counts", X = "SampleType"), "list")
    
    # Expect warning when over 10 taxa, expect error when over 20 taxa
    x_sub <- x[1:11]
    expect_warning(plotSeries(x_sub, abund_values = "counts", X = "SampleType"))
    x_sub <- x[1:21]
    expect_error(plotSeries(x_sub, abund_values = "counts", X = "SampleType"))
    
    ##################### Test .get_assay_data #################################
    x_sub <- x[1:10]
    expect_equal(miaViz:::.get_assay_data(x_sub, "counts"), assays(x_sub)$counts)
    
    ################### Test .incorporate_series_vis ###########################
    # Get data from colData
    X_test <- colData(x_sub)[, "SampleType"]
    colour_by_test <- colData(x_sub)[, "Primer"]
    X_and_colour_by <- data.frame(X = X_test, colour_by = colour_by_test)
    
    # Get data from rowData
    linetype_by_test <- rowData(x)[,"Family"]
    size_by_test <- rowData(x)[,"Kingdom"]
    linetype_by_and_size_by <- data.frame(linetype_by = linetype_by_test, size_by = size_by_test)
    
    # Get data from function
    data_from_function <- miaViz:::.incorporate_series_vis(x, X = "SampleType", 
                                                           colour_by = "Primer", linetype_by = "Family",
                                                           size_by = "Kingdom")
    
    # Divide list to 2 data frames
    sample_data_test <- data_from_function$sample_data
    feature_data_test <- data_from_function$feature_data
    
    # Expect that data frames are equal
    expect_equal(sample_data_test, X_and_colour_by)
    expect_equal(feature_data_test, linetype_by_and_size_by)
    
    ################## Test .melt_series_data ##################################
    melted <- miaViz:::.melt_series_data(assays(x)$counts, sample_data_test, feature_data_test, rownames(x))
    
    # Check 10 different random combinations
    for ( i in c(1:10) ){
        # Get random taxa name
        taxa <- rownames(x)[sample(length(rownames(x)), 1)]
        # Get random sample name
        sample <- colnames(x)[sample(length(colnames(x)), 1)]
        
        # Get value from assay data
        assay_value <- assays(x)$counts[taxa, sample]
        
        # Get value from melted data
        melted_assay_value <- melted[melted[, "feature"] == taxa,][melted[melted[, "feature"] == taxa,][, "sample"] == sample, ][, "Y"]
        
        ######### Expect that assay data is melted correctly ########
        expect_equal(melted_assay_value, assay_value)
        
        # Primer was chosen to be value of colour_by
        colData_colour_by_value <- colData(x)[sample, ][, "Primer"]
        melted_colour_by_value <- melted[melted[, "feature"] == taxa,][melted[melted[, "feature"] == taxa,][, "sample"] == sample, ][, "colour_by"]
        
        ######### Expect that colour_by values are melted correctly ########
        expect_equal(melted_colour_by_value, colData_colour_by_value)
        
        # Family was chosen to be value of linetype_by
        rowData_linetype_by_value <- rowData(x)[taxa, "Family"]
        melted_linetype_by_value <- melted[melted[, "feature"] == taxa,][melted[melted[, "feature"] == taxa,][, "sample"] == sample, ][, "linetype_by"]
        
        ######### Expect that linetype_by is melted correctly ########
        expect_equal(melted_linetype_by_value, rowData_linetype_by_value)
        
        # Kingdom was chosen to be value of size_by
        rowData_size_by_value <- rowData(x)[taxa, "Kingdom"]
        melted_size_by_value <- melted[melted[, "feature"] == taxa,][melted[melted[, "feature"] == taxa,][, "sample"] == sample, ][, "size_by"]
        
        ######### Expect that sizetype_by is melted correctly #########
        expect_equal(melted_size_by_value, rowData_size_by_value)
        
    }
    
})
