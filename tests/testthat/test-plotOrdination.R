
.make_tse <- function() {
    tse <- makeTSE(ncol = 2, nrow = 10)
    assayNames(tse) <- "counts"

    reducedDim(tse, "PCA") <- matrix(
        seq_len(6),
        ncol = 3,
        dimnames = list(colnames(tse), c("PC1", "PC2", "PC3"))
    )

    attr(reducedDim(tse, "PCA"), "eig") <- c(60, 30, 10)

    colData(tse)$numeric <- c(1, 2)
    colData(tse)$time <- c(1, 2)

    return(tse)
}

test_that("defaults are returned", {

    tse <- .make_tse()

    args <- .check_ordination_input(tse, "PCA")

    expect_identical(args$x, tse)
    expect_identical(args$dimred, "PCA")
    expect_identical(args$ncomponents, 1:2)
    expect_true(args$add.expl.var == FALSE)

})

test_that("scalar ncomponents expands", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        ncomponents = 2L
    )

    expect_identical(args$ncomponents, c(1L, 2L))

})

test_that("vector ncomponents is preserved", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        ncomponents = c(2L, 3L)
    )

    expect_identical(args$ncomponents, c(2L, 3L))

})

test_that("three components are rejected", {

    expect_error(
        .check_ordination_input(
            .make_tse(),
            "PCA",
            ncomponents = 1:3
        ),
        "ncomponents"
    )

})

test_that("component outside range is rejected", {

    expect_error(
        .check_ordination_input(
            .make_tse(),
            "PCA",
            ncomponents = c(1L, 4L)
        ),
        "ncomponents"
    )

})

test_that("non-integer components are rejected", {

    expect_error(
        .check_ordination_input(
            .make_tse(),
            "PCA",
            ncomponents = c(1, 2.5)
        ),
        "ncomponents"
    )

})

test_that("density and fill are mutually exclusive", {

    expect_error(
        .check_ordination_input(
            .make_tse(),
            "PCA",
            fill.by = "group",
            add.density = TRUE
        ),
        "Both 'add.density' and 'fill.by' cannot be specified simultaneously."
    )

})

for (arg in c(
    "add.points",
    "add.ellipse",
    "add.density",
    "add.centroids",
    "add.centroids.lines",
    "add.vectors",
    "add.rotation",
    "add.expl.var"
)) {

    test_that(paste(arg, "must be logical"), {

        x <- list(
            x = .make_tse(),
            dimred = "PCA"
        )

        x[[arg]] <- 1

        expect_error(
            do.call(.check_ordination_input, x),
            arg
        )

    })

}

test_that("extra arguments are propagated", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        point.alpha = .5,
        point.shape = 17
    )

    expect_identical(args$point.alpha, .5)
    expect_identical(args$point.shape, 17)

})

test_that("metadata arguments are propagated", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        colour.by = "group",
        fill.by = "group",
        shape.by = "group",
        size.by = "numeric",
        group.by = "group",
        linetype.by = "group",
        pair.by = "ID",
        sort.by = "time",
        facet.by = "group"
    )

    expect_identical(args$colour.by, "group")
    expect_identical(args$fill.by, "group")
    expect_identical(args$shape.by, "group")
    expect_identical(args$size.by, "numeric")
    expect_identical(args$group.by, "group")
    expect_identical(args$linetype.by, "group")
    expect_identical(args$pair.by, "ID")
    expect_identical(args$sort.by, "time")
    expect_identical(args$facet.by, "group")

})

test_that("plotOrdination returns a ggplot", {
    p <- plotOrdination(.make_tse(), "PCA")
    expect_s3_class(p, "ggplot")
})

test_that("integer ncomponents of length one expands to first two components", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        ncomponents = 2L
    )

    expect_identical(args$ncomponents, c(1L, 2L))

})

test_that("component order is preserved", {

    args <- .check_ordination_input(
        .make_tse(),
        "PCA",
        ncomponents = c(3L, 1L)
    )

    expect_identical(args$ncomponents, c(3L, 1L))

})

test_that("all plotting flags can be enabled", {

    expect_no_error(

        .check_ordination_input(
            .make_tse(),
            "PCA",
            add.points = TRUE,
            add.ellipse = TRUE,
            add.centroids = TRUE,
            add.centroids.lines = TRUE,
            add.vectors = TRUE,
            add.rotation = TRUE,
            add.expl.var = TRUE
        )

    )

})
