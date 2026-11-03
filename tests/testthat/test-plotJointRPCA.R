
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

test_that("vectors are hidden by default but can be added", {

    .n_segments <- function(p) sum(vapply(
        p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1)))

    expect_equal(.n_segments(plotJointRPCA(.make_joint(), "JointRPCA")), 0)

    expect_equal(.n_segments(plotJointRPCA(
        .make_joint(), "JointRPCA", add.vectors = TRUE)), 1)

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

############### plotJointRPCACovariance, ExplainedVariance, Network ############

# Joint-RPCA result with orthonormal loadings: 5 features in layer A and 4 in
# layer B, 3 components, 6 samples.
.make_joint_cov <- function() {
    set.seed(1)
    V <- qr.Q(qr(matrix(rnorm(9 * 3), ncol = 3)))
    dimnames(V) <- list(c(paste0("a", 1:5), paste0("b", 1:4)),
        paste0("PC", 1:3))
    rd <- matrix(rnorm(6 * 3), ncol = 3,
        dimnames = list(paste0("s", 1:6), paste0("PC", 1:3)))
    class(rd) <- c("JointRPCA", class(rd))
    attr(rd, "rotation") <- V
    attr(rd, "varExplained") <- c(PC1 = 3, PC2 = 2, PC3 = 1)
    attr(rd, "n_features") <- c(A = 5, B = 4)
    tse <- TreeSummarizedExperiment(assays = list(
        counts = matrix(0, 2, 6, dimnames = list(NULL, paste0("s", 1:6)))))
    reducedDim(tse, "JointRPCA") <- rd
    return(tse)
}

.make_joint_cov_mae <- function() {
    tse <- .make_joint_cov()
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(a = tse))
    metadata(mae)[["JointRPCA"]] <- reducedDim(tse, "JointRPCA")
    return(mae)
}

# Long data to matrix
.long_to_matrix <- function(df) {
    m <- xtabs(value ~ FeatureID + SampleID, df)
    unclass(m)[levels(df$FeatureID), levels(df$SampleID), drop = FALSE]
}

test_that("Joint-RPCA feature plots return ggplot for SCE and MAE", {
    for (x in list(.make_joint_cov(), .make_joint_cov_mae())) {
        expect_s3_class(plotJointRPCACovariance(x, "JointRPCA"), "ggplot")
        expect_s3_class(
            plotJointRPCACovariance(x, "JointRPCA", correlation = TRUE),
            "ggplot")
        expect_s3_class(
            plotJointRPCAExplainedVariance(x, "JointRPCA"), "ggplot")
        expect_s3_class(
            plotJointRPCAExplainedVariance(x, "JointRPCA", by = "feature"),
            "ggplot")
        expect_s3_class(plotJointRPCANetwork(x, "JointRPCA"), "ggplot")
    }
})

test_that("covariance equals V diag(s)^2 V^T", {
    res <- reducedDim(.make_joint_cov(), "JointRPCA")
    V <- attr(res, "rotation")
    s <- attr(res, "varExplained")
    expected <- V %*% diag(s^2) %*% t(V)
    df <- .get_jrpca_covariance_data(res, ntop = NULL, correlation = FALSE)
    m <- .long_to_matrix(df)
    expect_equal(unname(m), unname(expected[rownames(m), colnames(m)]))
    # Layers are added for rows and columns
    expect_equal(unname(df$Layer), ifelse(
        grepl("^a", df$FeatureID), "A", "B"))
    expect_equal(unname(df$ColLayer), ifelse(
        grepl("^a", df$SampleID), "A", "B"))
})

test_that("ntop selects features with the largest variance per layer", {
    res <- reducedDim(.make_joint_cov(), "JointRPCA")
    V <- attr(res, "rotation")
    s <- attr(res, "varExplained")
    variance <- rowSums(sweep(V^2, 2, s^2, "*"))
    keep <- .select_jrpca_features(res, ntop = 2)
    expect_length(keep, 4)
    expected <- c(
        names(sort(variance[1:5], decreasing = TRUE))[1:2],
        names(sort(variance[6:9], decreasing = TRUE))[1:2])
    expect_setequal(keep, expected)
    # ntop larger than layer keeps all features of the layer
    expect_length(.select_jrpca_features(res, ntop = 100), 9)
    expect_length(.select_jrpca_features(res, ntop = NULL), 9)
    # Same subset as calculating the full matrix and subsetting it
    full <- V %*% diag(s^2) %*% t(V)
    m <- .get_jrpca_matrix(res, ntop = 2, correlation = FALSE)
    expect_equal(m, full[keep, keep])
})

test_that("correlation is scaled covariance", {
    res <- reducedDim(.make_joint_cov(), "JointRPCA")
    m <- .get_jrpca_matrix(res, ntop = NULL, correlation = TRUE)
    expect_equal(unname(diag(m)), rep(1, 9))
    expect_true(all(abs(m) <= 1 + 1e-12))
    expect_equal(m, cov2cor(
        .get_jrpca_matrix(res, ntop = NULL, correlation = FALSE)))
})

test_that("explained variance sums to one per layer and feature", {
    res <- reducedDim(.make_joint_cov(), "JointRPCA")
    V <- attr(res, "rotation")
    s <- attr(res, "varExplained")
    # By layer
    df <- .get_jrpca_explained_variance_data(res, "layer", ntop = NULL)
    expect_setequal(as.character(df$FeatureID), c("A", "B"))
    expect_equal(nrow(df), 2 * 3)
    expect_equal(as.vector(tapply(df$value, df$FeatureID, sum)), c(1, 1))
    raw <- colSums(sweep(V[1:5, ]^2, 2, s^2, "*"))
    expect_equal(
        df$value[df$FeatureID == "A"][order(df$SampleID[df$FeatureID == "A"])],
        unname(raw / sum(raw)))
    # By feature
    df <- .get_jrpca_explained_variance_data(res, "feature", ntop = 2)
    expect_equal(nrow(df), 4 * 3)
    expect_equal(as.vector(tapply(df$value, df$FeatureID, sum)), rep(1, 4))
})

test_that("network keeps edges above threshold", {
    res <- reducedDim(.make_joint_cov(), "JointRPCA")
    m <- .get_jrpca_matrix(res, ntop = NULL, correlation = TRUE)
    data <- .get_jrpca_network_data(res, NULL, TRUE, threshold = 0.5)
    edges <- as.data.frame(tidygraph::activate(data$graph, "edges"))
    expect_true(all(abs(edges$value) >= 0.5))
    expect_equal(nrow(edges), sum(abs(m[upper.tri(m)]) >= 0.5))
    expect_equal(nrow(data$se), 9)
    expect_equal(unname(rowData(data$se)$Layer), rep(c("A", "B"), c(5, 4)))
    # No edges
    data <- .get_jrpca_network_data(res, NULL, TRUE, threshold = 2)
    expect_equal(igraph::ecount(data$graph), 0)
    expect_s3_class(
        plotJointRPCANetwork(.make_joint_cov(), "JointRPCA", threshold = 2),
        "ggplot")
    # Edge sign and strength are stored for plotting
    data <- .get_jrpca_network_data(res, NULL, TRUE, threshold = 0.5)
    edges <- as.data.frame(tidygraph::activate(data$graph, "edges"))
    expect_equal(edges$Sign,
        ifelse(edges$value > 0, "positive", "negative"))
    expect_equal(edges$Strength, abs(edges$value))
    # Nodes are shaped by layer, and the layers get different shapes
    p <- plotJointRPCANetwork(.make_joint_cov(), "JointRPCA",
        threshold = 0.5, show.label = TRUE)
    is_node <- vapply(
        p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
    expect_true("shape" %in% names(p$layers[[which(is_node)]]$mapping))
    node_data <- ggplot_build(p)$data[[which(is_node)]]
    expect_length(unique(node_data$shape), 2)
})

test_that("Joint-RPCA feature plots check input", {
    tse <- .make_joint_cov()
    expect_error(plotJointRPCACovariance(tse, "nope"), "dimred")
    expect_error(
        plotJointRPCACovariance(.make_joint_cov_mae(), "nope"), "dimred")
    expect_error(plotJointRPCACovariance(tse, "JointRPCA", ntop = 0), "ntop")
    expect_error(
        plotJointRPCACovariance(tse, "JointRPCA", ntop = "a"), "ntop")
    expect_error(
        plotJointRPCACovariance(tse, "JointRPCA", correlation = 1),
        "correlation")
    expect_error(
        plotJointRPCAExplainedVariance(tse, "JointRPCA", by = "sample"),
        "arg")
    expect_error(
        plotJointRPCANetwork(tse, "JointRPCA", threshold = "a"), "threshold")
    reducedDim(tse, "PCA") <- matrix(rnorm(12), ncol = 2)
    expect_error(plotJointRPCACovariance(tse, "PCA"), "JointRPCA")
    expect_error(plotJointRPCAExplainedVariance(tse, "PCA"), "JointRPCA")
    expect_error(plotJointRPCANetwork(tse, "PCA"), "JointRPCA")
})
