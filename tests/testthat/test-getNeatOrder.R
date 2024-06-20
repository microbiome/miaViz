testthat::set_max_fails(Inf)


context("neatsort")
neatsort_matrix <- matrix(c(10, 8,  2, 8,  3, 
                            4, 4,  4, 6,  5, 
                            5, 4,  4, 1,  2, 
                            2, 3,  5, 7,  8), nrow = 4, ncol = 5, byrow = TRUE, 
                          dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                          c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test getNeatOrder function", {
    # Test with valid inputs
    result <- getNeatOrder(neatsort_matrix[, c(1,2)], centering = "mean")
    expected <- c("Sample4", "Sample2", "Sample3", "Sample1")
    expect_equal(result, expected)
    
    # Test with no method input
    result <- getNeatOrder(neatsort_matrix[, c(1,2)])
    expected <- c("Sample4", "Sample2", "Sample3", "Sample1")
    expect_equal(result, expected)
})


context("check_args")
check_args_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2", "PC3", "PC4", "PC5")))

test_that("Test check_args method", {
    # Argument errors
    expect_error(.check_args(),
                 "argument \"x\" is missing, with no default")
    expect_error(.check_args(check_args_matrix[, c(1,2)]),
                 "argument \"centering\" is missing, with no default")
    
    # Non-matrix input
    expect_error(.check_args(as.data.frame(check_args_matrix)[, c(1,2)], centering = "mean"),
                 "Input data must be a matrix.")
    
    # Test for empty matrix
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    expect_error(.check_args(empty_matrix, centering = "mean"),
                 "No data to plot. Matrix must have at least one row and one column.")
    
    # Test for invalid number of columns
    expect_error(.check_args(check_args_matrix, centering = "mean"),
                 "Matrix must have only 2 columns.")
    
    # Test for invalid method
    expect_error(.check_args(check_args_matrix[, c(1,2)], "invalid_method"),
                 "'arg' should be one of “mean”, “median”, “none”")
    
    # Test for non-unique row names
    non_unique_row_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample1", "Sample4"),
                                                    c("PC1", "PC2", "PC3", "PC4", "PC5")))
    expect_error(.check_args(non_unique_row_matrix[, c(1,2)], centering = "mean"),
                 "Row names of the matrix must be unique.")
    
    # Test for non-unique column names
    non_unique_col_matrix <- matrix(1:20, nrow = 4, ncol = 5,
                                    dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                    c("PC1", "PC1", "PC3", "PC4", "PC5")))
    expect_error(.check_args(non_unique_col_matrix[, c(1,2)], centering = "mean"),
                 "Column names of the matrix must be unique.")
})


context("radial_theta")
radial_theta_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE,
                              dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                                             c("PC1", "PC2")))

test_that("Test radial_theta method", {
    # Argument errors
    expect_error(.radial_theta(),
                 "argument \"centering\" is missing, with no default")
    expect_error(.radial_theta(radial_theta_matrix),
                 "argument \"centering\" is missing, with no default")
    
    # Centering by mean
    centering_method <- "mean"
    result <- .radial_theta(radial_theta_matrix, centering_method)
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # Centering by median
    centering_method <- "median"
    result <- .radial_theta(radial_theta_matrix, centering_method)
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # No centering
    centering_method <- "none"
    result <- .radial_theta(radial_theta_matrix, centering_method)
    centered_data <- radial_theta_matrix
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
})


context("get_sorted_rownames")
theta_values <- c(Sample1 = 0.5, Sample2 = -1.2, Sample3 = 1.5, Sample4 = -0.8)

test_that("Test get_sorted_rownames method", {
    # Argument errors
    expect_error(.get_sorted_rownames(),
                 'argument "theta_values" is missing')
    
    # Valid sorting
    expected <- c("Sample2", "Sample4", "Sample1", "Sample3")
    result <- .get_sorted_rownames(theta_values)
    expect_equal(result, expected)
    
    # Edge case: all theta values are the same
    theta_values <- c(Sample1 = 0.5, Sample2 = 0.5, Sample3 = 0.5, Sample4 = 0.5)
    result <- .get_sorted_rownames(theta_values)
    expected <- c("Sample1", "Sample2", "Sample3", "Sample4")  # Order should be maintained
    expect_equal(result, expected)

    # Edge case: theta values contain NULL
    theta_values <- c(Sample1 = 0.5, Sample2 = NULL, Sample3 = 1.5, Sample4 = -0.8)
    result <- .get_sorted_rownames(theta_values)
    expected <- c("Sample4", "Sample1", "Sample3")
    expect_equal(result, expected)
})

