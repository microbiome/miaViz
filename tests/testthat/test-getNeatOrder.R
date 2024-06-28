context("neatsort")
neatsort_matrix <- matrix(c(10, 8,  2, 8,  3, 
                            4, 4,  4, 6,  5, 
                            5, 4,  4, 1,  2, 
                            2, 3,  5, 7,  8), nrow = 4, ncol = 5, byrow = TRUE, 
                          dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                          c("PC1", "PC2", "PC3", "PC4", "PC5")))

neatsort_matrix_indices <- matrix(c(10, 8,  2, 8,  3, 
                                    4, 4,  4, 6,  5, 
                                    5, 4,  4, 1,  2, 
                                    2, 3,  5, 7,  8), nrow = 4, ncol = 5, byrow = TRUE)

test_that("Test getNeatOrder function", {
    # Test with valid inputs
    result <- getNeatOrder(neatsort_matrix[, c(1,2)], centering = "mean")
    expected <- c(4, 2, 3, 1)
    expect_equal(result, expected)
    
    # Test with no method input
    result <- getNeatOrder(neatsort_matrix[, c(1,2)])
    expected <- c(4, 2, 3, 1)
    expect_equal(result, expected)
    
    # Test with indice matrix
    result <- getNeatOrder(neatsort_matrix_indices[, c(1,2)], centering = "mean")
    expected <- c(4, 2, 3, 1)
    expect_equal(result, expected_indices)
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
    error_message <- expect_error(.check_args(check_args_matrix[, c(1,2)], "invalid_method"))
    expect_true(grepl("'arg' should be one of", error_message))
})


context("radial_theta")
radial_theta_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE,
                              dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                                                             c("PC1", "PC2")))

radial_theta_matrix_indices <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE)

test_that("Test radial_theta method", {
    # Argument errors
    expect_error(.radial_theta(),
                 "argument \"centering\" is missing, with no default")
    expect_error(.radial_theta(radial_theta_matrix),
                 "argument \"centering\" is missing, with no default")
    
    # Centering by mean
    result <- .radial_theta(radial_theta_matrix, "mean")
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # Centering by median
    result <- .radial_theta(radial_theta_matrix, "median")
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE,
                            dimnames = list(c("Sample1", "Sample2", "Sample3", "Sample4"),
                                            c("PC1", "PC2")))
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # No centering
    result <- .radial_theta(radial_theta_matrix, "none")
    centered_data <- radial_theta_matrix
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    names(expected) <- rownames(centered_data)
    expect_equal(result, expected)
    
    # Test with indice matrix
    result <- .radial_theta(radial_theta_matrix_indices, "mean")
    centered_data <- matrix(c(-3, -3, -1, -1, 1, 1, 3, 3), ncol = 2, byrow = TRUE)
    expected <- atan2(centered_data[, 2], centered_data[, 1])
    expect_equal(result, expected)
})
