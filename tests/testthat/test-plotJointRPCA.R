
.make_joint <- function() {

    tse <- makeTSE(nrow = 8, ncol = 4)
    assayNames(tse) <- "counts"

    rd <- matrix(rnorm(8), ncol = 2)
    class(rd) <- c("JointRPCA", class(rd))

    rotation <- matrix(rnorm(16), ncol = 2)
    rownames(rotation) <- paste0("feature", seq_len(nrow(rotation)))

    attr(rd, "rotation") <- rotation
    attr(rd, "n_features") <- c(Layer1 = 4, Layer2 = 4)

    reducedDim(tse, "JointRPCA") <- rd

    return(tse)
}

test_that("returns ggplot", {

    p <- plotJointRPCA(.make_joint(), "JointRPCA")

    expect_s3_class(p, "ggplot")

})

test_that("FALSE disables vectors", {

    expect_no_error(

        plotJointRPCA(
            .make_joint(),
            "JointRPCA",
            add.vectors = FALSE
        )

    )

})

test_that("character vectors accepted", {

    expect_no_error(

        plotJointRPCA(
            .make_joint(),
            "JointRPCA",
            add.vectors = "feature"
        )

    )

})

test_that("invalid add.vectors errors", {

    expect_error(

        plotJointRPCA(
            .make_joint(),
            "JointRPCA",
            add.vectors = 1
        ),

        "add.vectors"

    )

})

test_that("non JointRPCA reducedDim errors", {

    tse <- makeTSE(nrow = 5, ncol = 2)

    reducedDim(tse,"PCA") <- matrix(rnorm(4),2)

    expect_error(

        plotJointRPCA(
            tse,
            "PCA"
        ),

        "JointRPCA"

    )

})

test_that("returns data.frame", {

    x <- .get_joint_rpca_vector_data(
        .make_joint(),
        "JointRPCA",
        TRUE
    )

    expect_s3_class(x,"data.frame")

})

test_that("FALSE returns NULL", {

    expect_null(

        .get_joint_rpca_vector_data(
            .make_joint(),
            "JointRPCA",
            FALSE
        )

    )

})

test_that("feature filtering works", {

    x <- .get_joint_rpca_vector_data(
        .make_joint(),
        "JointRPCA",
        "feature1"
    )

    expect_true(

        all(grepl(
            "feature1",
            x$vector_label
        ))

    )

})

test_that("case insensitive filtering works", {

    x <- .get_joint_rpca_vector_data(
        .make_joint(),
        "JointRPCA",
        "FEATURE",
        ignore.case = TRUE
    )

    expect_gt(nrow(x),0)

})

test_that("ntop limits each layer", {

    x <- .get_joint_rpca_vector_data(
        .make_joint(),
        "JointRPCA",
        TRUE,
        ntop = 2
    )

    expect_true(

        all(table(x$Layer) <= 2)

    )

})

test_that("ntop NULL keeps all", {

    x <- .get_joint_rpca_vector_data(
        .make_joint(),
        "JointRPCA",
        TRUE,
        ntop = NULL
    )

    expect_equal(

        nrow(x),

        8

    )

})

test_that("ignore.case validated", {

    expect_error(

        .get_joint_rpca_vector_data(
            .make_joint(),
            "JointRPCA",
            TRUE,
            ignore.case = 1
        ),

        "ignore.case"

    )

})

test_that("ntop validated", {

    expect_error(

        .get_joint_rpca_vector_data(
            .make_joint(),
            "JointRPCA",
            TRUE,
            ntop = "a"
        ),

        "ntop"

    )

})

.make_vectors <- function() {

    data.frame(
        PC1 = c(.2,.5),
        PC2 = c(.4,.1),
        Layer = c("A","B"),
        vector_label = c("x","y")
    )

}

test_that("returns ggplot", {

    p <- ggplot()

    p2 <- .add_ordination_vectors(
        p,
        .make_vectors()
    )

    expect_s3_class(p2,"ggplot")

    expect_length(

        p2$layers,

        2

    )

})

test_that("supports geom_label", {

    p <- .add_ordination_vectors(
        ggplot(),
        .make_vectors(),
        text.labels = FALSE
    )

    expect_s3_class(p,"ggplot")

})

test_that("supports non-repel labels", {

    p <- .add_ordination_vectors(
        ggplot(),
        .make_vectors(),
        repel.labels = FALSE
    )

    expect_s3_class(p,"ggplot")

})

test_that("supports label without repel", {

    p <- .add_ordination_vectors(
        ggplot(),
        .make_vectors(),
        repel.labels = FALSE,
        text.labels = FALSE
    )

    expect_s3_class(p,"ggplot")

})
