testthat::set_max_fails(Inf)


context("take_subset")
take_subset_matrix <- matrix(1:20, ncol = 5, byrow = TRUE, dimnames = list(c("Feature1", "Feature2", "Feature3", "Feature4"), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")))

test_that("Test take_subset method", {
    # Argument errors
    expect_error(miaViz:::.take_subset(),
                 'argument "data" is missing')
    
    # Valid subsetting by row names
    subset <- c("Feature1", "Feature3", "Feature4")
    result <- miaViz:::.take_subset(take_subset_matrix, subset)
    expected <- matrix(c(1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20), nrow = 3, byrow = TRUE,
                       dimnames = list(c("Feature1", "Feature3", "Feature4"), 
                                       c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")))
    expect_equal(result, expected)
    
    # Valid subsetting by row indices
    subset <- c(1, 3, 4)
    result <- miaViz:::.take_subset(take_subset_matrix, subset)
    expected <- matrix(c(1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20), nrow = 3, byrow = TRUE,
                       dimnames = list(c("Feature1", "Feature3", "Feature4"), 
                                       c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")))
    expect_equal(result, expected)
    
    # Subsetting with indices out of bounds
    subset <- c(1, 6)
    expect_error(miaViz:::.take_subset(take_subset_matrix, subset), "subscript out of bounds")
    
    # Subsetting with non-existent row names
    subset <- c("Feature1", "Feature6")
    expect_error(miaViz:::.take_subset(take_subset_matrix, subset), "subscript out of bounds")
    
    # Subsetting with incorrect types
    subset <- list(1, 3)
    expect_error(miaViz:::.take_subset(take_subset_matrix, subset), "invalid subscript type")
})


context("take_dimensions")
take_dimensions_matrix <- matrix(1:20, ncol = 5, byrow = TRUE, dimnames = list(c("Feature1", "Feature2", "Feature3", "Feature4"), c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test take_dimensions method",{
    # Argument errors
    expect_error(miaViz:::.take_dimensions(),
                 'argument "data" is missing')
    
    # Valid subsetting by column indices
    dimensions <- c(1, 3)
    result <- miaViz:::.take_dimensions(take_dimensions_matrix, dimensions)
    expected <- matrix(c(1, 6, 11, 16, 3, 8, 13, 18), nrow = 4, byrow = FALSE,
                       dimnames = list(c("Feature1", "Feature2", "Feature3", "Feature4"), 
                                       c("PC1", "PC3")))
    expect_equal(result, expected)
    
    # Valid subsetting by column names
    dimensions <- c("PC1", "PC3")
    result <- miaViz:::.take_dimensions(take_dimensions_matrix, dimensions)
    expected <- matrix(c(1, 6, 11, 16, 3, 8, 13, 18), nrow = 4, byrow = FALSE,
                       dimnames = list(c("Feature1", "Feature2", "Feature3", "Feature4"), 
                                       c("PC1", "PC3")))
    expect_equal(result, expected)
    
    # Subsetting with indices out of bounds
    subset <- c(1, 6)
    expect_error(miaViz:::.take_dimensions(take_dimensions_matrix, subset), "subscript out of bounds")
    
    # Subsetting with non-existent row names
    dimensions <- c("PC1", "PC6")
    expect_error(miaViz:::.take_dimensions(take_dimensions_matrix, subset), "subscript out of bounds")
    
    # Subsetting with incorrect types
    subset <- list(1, 3)
    expect_error(miaViz:::.take_dimensions(take_dimensions_matrix, subset), "invalid subscript type")
})


context("radial_theta")
radial_theta_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE,
                              dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                                             c("PC1", "PC2")))

test_that("Test radial_theta method", {
    # Argument errors
    expect_error(miaViz:::.take_dimensions(),
                 'argument "data" is missing')
    
    # Centering by mean
    centering_method <- "mean"
    result <- miaViz:::.radial_theta(radial_theta_matrix, centering_method)
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # Centering by median
    centering_method <- "median"
    result <- miaViz:::.radial_theta(radial_theta_matrix, centering_method)
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # No centering
    centering_method <- "none"
    result <- miaViz:::.radial_theta(radial_theta_matrix, centering_method)
    centered_data <- radial_theta_matrix
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # Unsupported centering method
    centering_method <- "unsupported"
    expect_error(miaViz:::.radial_theta(radial_theta_matrix, centering_method), "Unsupported centering method. Choose either 'mean', 'median', 'mode', or 'none'.")
})


context("get_sorted_rownames")
theta_values <- c(Sample1 = 0.5, Sample2 = -1.2, Sample3 = 1.5, Sample4 = -0.8)

test_that("Test get_sorted_rownames method", {
    # Argument errors
    expect_error(miaViz:::.get_sorted_rownames(),
                 'argument "theta_values" is missing')
    expect_error(miaViz:::.get_sorted_rownames(theta_values),
                 'argument "sorting_order" is missing')
    
    # Valid sorting in ascending order
    sorting_order <- "ascending"
    expected <- c("Sample2", "Sample4", "Sample1", "Sample3")
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expect_equal(result, expected)
    
    # Sorting in descending order
    sorting_order <- "descending"
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expected <- c("Sample3", "Sample1", "Sample4", "Sample2")
    expect_equal(result, expected)
    
    # Edge case: all theta values are the same
    theta_values <- c(Sample1 = 0.5, Sample2 = 0.5, Sample3 = 0.5, Sample4 = 0.5)
    sorting_order <- "ascending"
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expected <- c("Sample1", "Sample2", "Sample3", "Sample4")  # Order should be maintained
    expect_equal(result, expected)
    
    sorting_order <- "descending"
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expected <- c("Sample1", "Sample2", "Sample3", "Sample4")  # Order should be maintained
    expect_equal(result, expected)
    
    # Edge case: theta values contain NULL
    theta_values <- c(Sample1 = 0.5, Sample2 = NULL, Sample3 = 1.5, Sample4 = -0.8)
    sorting_order <- "ascending"
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expected <- c("Sample4", "Sample1", "Sample3")
    expect_equal(result, expected)
    
    sorting_order <- "descending"
    result <- miaViz:::.get_sorted_rownames(theta_values, sorting_order)
    expected <- c("Sample3", "Sample1", "Sample4") 
    expect_equal(result, expected)
})


context("check_neatsort_args")
check_args_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                        dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                        c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test check_neatsort_args method", {
    # Argument errors
    expect_error(miaViz:::check_neatsort_args(),
                 "object 'check_neatsort_args' not found")
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix),
                 'argument "subset" is missing')
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, c(1, 3, 4)),
                 'argument \"dimensions\" is missing, with no default')
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, c(1, 3, 4), c(1, 3)),
                 'argument \"centering_method\" is missing, with no default')
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, c(1, 3, 4), c(1, 3), "mean"),
                 'argument "sorting_order" is missing')
    
    # Non-matrix input
    expect_error(miaViz:::.check_neatsort_args(as.data.frame(check_args_matrix), NULL, c(1, 2), "mean", "ascending"),
                 "Input data must be a matrix.")
    
    # Test for empty matrix
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    expect_error(miaViz:::.check_neatsort_args(empty_matrix, NULL, c(1, 2), "mean", "ascending"),
                 "No data to plot. Matrix must have at least one row and one column.")
    
    # Test valid subset
    expect_silent(miaViz:::.check_neatsort_args(check_args_matrix, NULL, c(1, 2), "mean", "ascending"))
    expect_silent(miaViz:::.check_neatsort_args(check_args_matrix, c(1, 3, 4), c(1, 2), "mean", "ascending"))
    expect_silent(miaViz:::.check_neatsort_args(check_args_matrix, c("Sample1", "Sample2"), c(1, 2), "mean", "ascending"))
    
    # Test invalid subset
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, 5, c(1, 2), "mean", "ascending"),
                 "Subset refers to rows that do not exist in the data.")
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, "Sample5", c(1, 2), "mean", "ascending"),
                 "Subset refers to row names that do not exist in the data.")
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, list(1, 2), c(1, 2), "mean", "ascending"),
                 "Subset must be a vector of row indices or names.")
    
    # Test for invalid dimensions
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, NULL, c(1, 6), "mean", "ascending"),
                 "dimensions refer to columns that do not exist in the data.")
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, NULL, c(1, 2, 3), "mean", "ascending"),
                 "Exactly two dimensions must be specified.")
    
    # Test for invalid centering_method
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, NULL, c(1, 2), "invalid_method", "ascending"),
                 "'arg' should be one of “mean”, “median”, “none”")
    
    # Test for invalid sorting_order
    expect_error(miaViz:::.check_neatsort_args(check_args_matrix, NULL, c(1, 2), "mean", "invalid_order"),
                 "'arg' should be one of “ascending”, “descending”")
    
    # Test for non-unique row names
    non_unique_row_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample1", "Sample4"),
                                                    c("PC1", "PC2", "PC3", "PC4", "PC5")))
    expect_error(miaViz:::.check_neatsort_args(non_unique_row_matrix, NULL, c(1, 2), "mean", "ascending"),
                 "Row names of the matrix must be unique.")
    
    # Test for non-unique column names
    non_unique_col_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                    c("PC1", "PC2", "PC1", "PC4", "PC5")))
    expect_error(miaViz:::.check_neatsort_args(non_unique_col_matrix, NULL, c(1, 2), "mean", "ascending"),
                 "Column names of the matrix must be unique.")
    
})


context("neatsort")
neatsort_matrix <- matrix(c(10, 8,  2, 8,  3, 
                            4, 4,  4, 6,  5, 
                            5, 4,  4, 1,  2, 
                            2, 3,  5, 7,  8), nrow = 4, ncol = 5, byrow = TRUE, 
                          dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                          c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test neatsort function", {
    # Test with valid inputs, no subset
    result <- miaViz::neatsort(neatsort_matrix, subset = NULL, dimensions = c(1, 2), centering_method = "mean", sorting_order = "ascending")
    expected <- c("Sample4", "Sample2", "Sample3", "Sample1")
    expect_equal(result, expected)
    
    # Test with valid inputs and subset
    subset <- c("Sample1", "Sample3")
    result <- miaViz::neatsort(neatsort_matrix, subset = subset, dimensions = c(1, 2), centering_method = "mean", sorting_order = "ascending")
    expected <- c("Sample3", "Sample1")
    expect_equal(result, expected)
})

