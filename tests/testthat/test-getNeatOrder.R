testthat::set_max_fails(Inf)


context("neatsort")
neatsort_matrix <- matrix(c(10, 8,  2, 8,  3, 
                            4, 4,  4, 6,  5, 
                            5, 4,  4, 1,  2, 
                            2, 3,  5, 7,  8), nrow = 4, ncol = 5, byrow = TRUE, 
                          dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                          c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test getNeatOrder function", {
    # Test with valid inputs, no subset
    result <- miaViz::getNeatOrder(neatsort_matrix, dimensions = c(1, 2), centering_method = "mean")
    expected <- c("Sample4", "Sample2", "Sample3", "Sample1")
    expect_equal(result, expected)
    
    # Test with valid inputs, no subset, decreasing parameter = TRUE
    result <- miaViz::getNeatOrder(neatsort_matrix, dimensions = c(1, 2), centering_method = "mean", decreasing = TRUE)
    expected <- c("Sample1", "Sample3", "Sample2", "Sample4")
    expect_equal(result, expected)
    
    # Test with dimensions col names
    result <- miaViz::getNeatOrder(neatsort_matrix, dimensions = c("PC1", "PC2"), centering_method = "mean")
    expected <- c("Sample4", "Sample2", "Sample3", "Sample1")
    expect_equal(result, expected)
    
    # Test with valid inputs and subset row names
    subset <- c("Sample1", "Sample3")
    result <- miaViz::getNeatOrder(neatsort_matrix, subset, dimensions = c(1, 2), centering_method = "mean")
    expected <- c("Sample3", "Sample1")
    expect_equal(result, expected)
    
    # Test with valid inputs and subset row indices
    subset <- c(1, 3)
    result <- miaViz::getNeatOrder(neatsort_matrix, subset, dimensions = c(1, 2), centering_method = "mean")
    expected <- c("Sample3", "Sample1")
    expect_equal(result, expected)
})


context("check_args")
check_args_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test check_args method", {
    # Argument errors
    expect_error(miaViz:::check_args(),
                 "object 'check_args' not found")
    expect_error(miaViz:::.check_args(check_args_matrix),
                 'argument "subset" is missing')
    expect_error(miaViz:::.check_args(check_args_matrix, c(1, 3, 4)),
                 'argument \"dimensions\" is missing, with no default')
    expect_error(miaViz:::.check_args(check_args_matrix, c(1, 3, 4), c(1, 3)),
                 'argument \"centering_method\" is missing, with no default')
    expect_error(miaViz:::.check_args(check_args_matrix, c(1, 3, 4), c(1, 3), "mean"),
                 'argument "decreasing" is missing')
    
    # Non-matrix input
    expect_error(miaViz:::.check_args(as.data.frame(check_args_matrix), NULL, c(1, 2), "mean", FALSE),
                 "Input data must be a matrix.")
    
    # Test for empty matrix
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    expect_error(miaViz:::.check_args(empty_matrix, NULL, c(1, 2), "mean", FALSE),
                 "No data to plot. Matrix must have at least one row and one column.")
    
    # Test valid subset
    expect_silent(miaViz:::.check_args(check_args_matrix, NULL, c(1, 2), "mean", FALSE))
    expect_silent(miaViz:::.check_args(check_args_matrix, c(1, 3, 4), c(1, 2), "mean", FALSE))
    expect_silent(miaViz:::.check_args(check_args_matrix, c("Sample1", "Sample2"), c(1, 2), "mean", FALSE))
    
    # Test invalid subset
    expect_error(miaViz:::.check_args(check_args_matrix, 5, c(1, 2), "mean", FALSE),
                 "Subset refers to rows that do not exist in the data.")
    expect_error(miaViz:::.check_args(check_args_matrix, "Sample5", c(1, 2), "mean", FALSE),
                 "Subset refers to row names that do not exist in the data.")
    expect_error(miaViz:::.check_args(check_args_matrix, list(1, 2), c(1, 2), "mean", FALSE),
                 "Subset must be a vector of row indices or names.")
    
    # Test for invalid dimensions
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, c(1, 6), "mean", FALSE),
                 "dimensions refer to columns that do not exist in the data.")
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, c("PC6"), "mean", FALSE),
                 "dimensions refer to column names that do not exist in the data.")
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, list("PC6"), "mean", FALSE),
                 "dimensions must be a vector of column indices or names.")
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, c(1, 2, 3), "mean", FALSE),
                 "Exactly two dimensions must be specified.")
    
    # Test for invalid centering_method
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, c(1, 2), "invalid_method", FALSE),
                 "'arg' should be one of “mean”, “median”, “none”")
    
    # Test for invalid sorting_order
    expect_error(miaViz:::.check_args(check_args_matrix, NULL, c(1, 2), "mean", "invalid_order"),
                 "decreasing must be a single boolean value.")
    
    # Test for non-unique row names
    non_unique_row_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample1", "Sample4"),
                                                    c("PC1", "PC2", "PC3", "PC4", "PC5")))
    expect_error(miaViz:::.check_args(non_unique_row_matrix, NULL, c(1, 2), "mean", FALSE),
                 "Row names of the matrix must be unique.")
    
    # Test for non-unique column names
    non_unique_col_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                    c("PC1", "PC2", "PC1", "PC4", "PC5")))
    expect_error(miaViz:::.check_args(non_unique_col_matrix, NULL, c(1, 2), "mean", FALSE),
                 "Column names of the matrix must be unique.")
    
})


context("radial_theta")
radial_theta_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE,
                              dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                                             c("PC1", "PC2")))

test_that("Test radial_theta method", {
    # Argument errors
    expect_error(miaViz:::.radial_theta(),
                 'argument \"centering_method\" is missing, with no default')
    expect_error(miaViz:::.radial_theta(radial_theta_matrix),
                 'argument "centering_method" is missing')
    
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
                 'argument "decreasing" is missing')
    
    # Valid sorting in ascending order
    expected <- c("Sample2", "Sample4", "Sample1", "Sample3")
    result <- miaViz:::.get_sorted_rownames(theta_values, FALSE)
    expect_equal(result, expected)
    
    # Sorting in descending order
    result <- miaViz:::.get_sorted_rownames(theta_values, TRUE)
    expected <- c("Sample3", "Sample1", "Sample4", "Sample2")
    expect_equal(result, expected)
    
    # Edge case: all theta values are the same
    theta_values <- c(Sample1 = 0.5, Sample2 = 0.5, Sample3 = 0.5, Sample4 = 0.5)
    result <- miaViz:::.get_sorted_rownames(theta_values, FALSE)
    expected <- c("Sample1", "Sample2", "Sample3", "Sample4")  # Order should be maintained
    expect_equal(result, expected)
    
    result <- miaViz:::.get_sorted_rownames(theta_values, TRUE)
    expected <- c("Sample1", "Sample2", "Sample3", "Sample4")  # Order should be maintained
    expect_equal(result, expected)
    
    # Edge case: theta values contain NULL
    theta_values <- c(Sample1 = 0.5, Sample2 = NULL, Sample3 = 1.5, Sample4 = -0.8)
    result <- miaViz:::.get_sorted_rownames(theta_values, FALSE)
    expected <- c("Sample4", "Sample1", "Sample3")
    expect_equal(result, expected)
    
    result <- miaViz:::.get_sorted_rownames(theta_values, TRUE)
    expected <- c("Sample3", "Sample1", "Sample4") 
    expect_equal(result, expected)
})

