#' @name
#' plotJointRPCA
#'
#' @title
#' Visualize Joint-RPCA results
#'
#' @description
#' These functions visualize Joint-RPCA results stored in a
#' \code{SingleCellExperiment} or \code{MultiAssayExperiment}.
#'
#' \itemize{
#'   \item \code{plotJointRPCA()} creates a two-dimensional ordination plot
#'   of the samples.
#'   \item \code{plotJointRPCACovariance()} creates a heatmap of the
#'   feature-feature covariance (or correlation), split into blocks by data
#'   layer.
#'   \item \code{plotJointRPCAExplainedVariance()} creates a heatmap showing
#'   how the variance of each data layer, or of each feature, is divided
#'   across the components.
#'   \item \code{plotJointRPCANetwork()} creates a network where nodes are
#'   features and edges connect features with strong covariance (or
#'   correlation).
#' }
#'
#' @details
#' \code{plotJointRPCA()} is a wrapper around
#' \code{\link[=plotOrdination]{plotOrdination}} for Joint-RPCA results.
#' Consequently, most graphical parameters are passed directly to
#' \code{plotOrdination()}, including options for colouring, grouping,
#' faceting and adding ellipses. See
#' \code{\link[=plotOrdination]{plotOrdination}} for a complete description of
#' these arguments.
#'
#' The feature covariance follows the gemelli implementation of Joint-RPCA:
#' \deqn{\Sigma = V \textrm{diag}(s)^2 V^T,}
#' where \eqn{V} are the joint feature loadings
#' (\code{attr(res, "rotation")}) and \eqn{s} the singular values
#' (\code{attr(res, "varExplained")}). The full matrix has one row and column
#' per feature and can be very large. Therefore, it is calculated only for the
#' \code{ntop} features with the largest variance (the diagonal of the
#' covariance matrix) from each layer. Features are ordered by hierarchical
#' clustering. The covariance is estimated from a low-rank representation.
#' With a small number of components, correlations are often close to -1 or 1,
#' and they reflect how features align in the Joint-RPCA space rather than
#' correlations in the observed data. The network uses the same covariance and
#' is plotted with \code{\link[=plotGraph]{plotRowGraph}}.
#'
#' In \code{plotJointRPCAExplainedVariance()}, the variance of feature
#' \eqn{i} on component \eqn{j} is \eqn{s_j^2 V_{ij}^2}. For a layer, the
#' values are summed over its features. Each row (layer or feature) is scaled
#' to sum to one. The values describe the low-rank Joint-RPCA model, not the
#' proportion of variance explained in the observed data.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}}
#' object containing Joint-RPCA results.
#'
#' @param dimred \code{character(1)}: Specifies the name of the Joint-RPCA
#' result to plot.
#'
#' @param ntop \code{integer(1)} or \code{NULL}. Number of features shown from
#' each data layer. The features with the largest variance are retained. Set
#' to \code{NULL} to show all features. In \code{plotJointRPCAExplainedVariance}
#' used only when \code{by = "feature"}.
#' (Default: \code{20} for \code{plotJointRPCACovariance} and
#' \code{plotJointRPCAExplainedVariance}, \code{10} for
#' \code{plotJointRPCANetwork})
#'
#' @param correlation \code{logical(1)}. Should correlations be shown instead
#' of covariances? (Default: \code{FALSE} for \code{plotJointRPCACovariance},
#' \code{TRUE} for \code{plotJointRPCANetwork})
#'
#' @param by \code{character(1)}. Whether to show variance by data layer
#' (\code{"layer"}) or by feature (\code{"feature"}).
#' (Default: \code{"layer"})
#'
#' @param threshold \code{numeric(1)}. Network edges with absolute value below
#' the threshold are dropped. (Default: \code{0.9})
#'
#' @param ...
#' Additional arguments.
#' \itemize{
#' \item \code{plotJointRPCA}: passed to
#' \code{\link[=plotOrdination]{plotOrdination}}.
#' \item \code{plotJointRPCANetwork}: passed to
#' \code{\link[=plotGraph]{plotRowGraph}}, for example \code{layout} and
#' \code{show.label}.
#' \item \code{plotJointRPCACovariance} and
#' \code{plotJointRPCAExplainedVariance}: currently not used.
#' }
#'
#' @examples
#' data("HintikkaXOData")
#' mae <- HintikkaXOData
#'
#' mae[[1]] <- transformAssay(
#'     mae[[1]],
#'     assay.type = "counts",
#'     method = "rclr",
#'     impute = FALSE
#' )
#' mae[[2]] <- transformAssay(
#'     mae[[2]],
#'     assay.type = "nmr",
#'     method = "log10"
#' )
#'
#' mae <- addJointRPCA(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("rclr", "log10")
#' )
#'
#' # Basic Joint-RPCA plot
#' plotJointRPCA(mae, "JointRPCA")
#'
#' # Colour samples by metadata
#' plotJointRPCA(mae, "JointRPCA", colour.by = "Fat")
#'
#' # Add confidence ellipses
#' plotJointRPCA(
#'     mae,
#'     "JointRPCA",
#'     colour.by = "Fat",
#'     add.ellipse = TRUE
#' )
#'
#' # Covariance of the top 20 features from each layer
#' plotJointRPCACovariance(mae, "JointRPCA")
#'
#' # Correlation of the top 10 features from each layer
#' plotJointRPCACovariance(mae, "JointRPCA", ntop = 10, correlation = TRUE)
#'
#' # Variance explained by each component in each layer
#' plotJointRPCAExplainedVariance(mae, "JointRPCA")
#'
#' # Variance explained by each component in each feature
#' plotJointRPCAExplainedVariance(mae, "JointRPCA", by = "feature")
#'
#' # Network of the top 10 features from each layer
#' plotJointRPCANetwork(mae, "JointRPCA")
#'
#' # Stricter threshold and node labels
#' plotJointRPCANetwork(mae, "JointRPCA", threshold = 0.97, show.label = TRUE)
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[=plotOrdination]{plotOrdination}}
#'   \item \code{\link[=plotGraph]{plotRowGraph}}
#' }
#'
NULL

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCA", signature = c(x = "MultiAssayExperiment"),
    function(x, dimred, ...){
        .check_metadata_present(dimred, x)
        # Construct TreeSE from results so that we can use TreeSE function
        # (and thus plotOrdination)
        res <- metadata(x)[[dimred]]
        #
        col_data <- x |> colData()
        col_data <- col_data[
            match(rownames(res), rownames(col_data)), , drop = FALSE]
        rownames(col_data) <- rownames(res)
        #
        tse <- TreeSummarizedExperiment(
            assays = SimpleList(counts = matrix(
                ncol = nrow(res), dimnames = list(NULL, rownames(res)))),
            colData = col_data,
            reducedDims = list(JointRPCA = res)
        )
        # Plot
        p <- plotJointRPCA(tse, dimred = "JointRPCA", ...)
        return(p)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCA", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ...){
        .check_dimred_present(dimred, x)
        if( !inherits(reducedDim(x, dimred), "JointRPCA") ){
            stop("The plotted results must be a class 'JointRPCA'.",
                call. = FALSE)
        }
        p <- .plot_joint_rpca_ordination(x, dimred, ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# Plot ordination and optionally feature loading vectors. Vectors are
# supported but not documented, as they do not look good enough yet. Enable
# with add.vectors = TRUE (or a character vector of feature name patterns).
# add.vectors is captured here so that it is not passed to plotOrdination(),
# where it has a different meaning.
.plot_joint_rpca_ordination <- function(x, dimred, add.vectors = FALSE, ...){
    if( !( .is_a_bool(add.vectors) || is.character(add.vectors)) ){
        stop("'add.vectors must be TRUE or FALSE or character vector.",
            call. = FALSE)
    }
    p <- plotOrdination(x, dimred, ...)
    vector_data <- .get_joint_rpca_vector_data(x, dimred, add.vectors, ...)
    if( !is.null(vector_data) ){
        p <- .add_ordination_vectors(p, vector_data, ...)
    }
    return(p)
}

# Get feature loadings that will be plotted as vectors
#' @importFrom dplyr group_by slice_max ungroup
.get_joint_rpca_vector_data <- function(
        tse, reduced_dim, add.vectors, ignore.case = FALSE, ntop = 10, ...){
    if( !.is_a_bool(ignore.case) ){
        stop("'ignore.case' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !(.is_an_integer(ntop) || is.null(ntop)) ){
        stop("'ntop' must be an integer.", call. = FALSE)
    }
    #
    # Get vector data, i.e, rotations or loadings
    res <- reducedDim(tse, reduced_dim)
    vector_data <- if(!(.is_a_bool(add.vectors) && !add.vectors))
        attributes(res)[["rotation"]]
    # If user wanted to plot them, wrangle the data
    if( !is.null(vector_data) ){
        # Add labels and layer name
        vector_data <- as.data.frame(vector_data)
        vector_data[["vector_label"]] <- rownames(vector_data)
        n_features <- attr(res, "n_features")
        vector_data[["Layer"]] <- rep(
            names(n_features),
            times = unname(n_features)
        )

        # Subset vectors by selecting only those ones that user has specified
        if( is.character(add.vectors) ){
            add.vectors <- paste0(add.vectors, collapse = "|")
            keep <- vapply(rownames(vector_data), function(x)
                grepl(add.vectors, x, perl = TRUE, ignore.case = ignore.case),
                logical(1L))
            vector_data <- vector_data[keep, ]
        }

        # Keep only the n longest vectors per layer. The idea is to show only
        # features that are highly associated with the showed ordination space.
        if (!is.null(ntop) && nrow(vector_data) > 0L) {
            xvar <- colnames(vector_data)[1]
            yvar <- colnames(vector_data)[2]

            vector_data[["length_temporary"]] <-
                sqrt(vector_data[[xvar]]^2 + vector_data[[yvar]]^2)

            vector_data <- vector_data |>
                group_by(Layer) |>
                slice_max(length_temporary, n = ntop, with_ties = FALSE) |>
                ungroup()

            vector_data[["length_temporary"]] <- NULL
        }

        # If all vectors were removed, give NULL
        if( nrow(vector_data) == 0L ){
            vector_data <- NULL
        }
    }

    return(vector_data)
}

# Add feature loadings as vectors
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom grid arrow unit
.add_ordination_vectors <- function(
        p,
        vector_data,
        vec.size = 0.5,
        arrow.size = 0.25,
        label.size = 4,
        vec.color = "black",
        label.color = "black",
        vec.text = TRUE,
        text.labels = TRUE,
        repel.labels = TRUE,
        min.segment.length = 0.5,
        box.padding = 0.25,
        point.padding = 1e-06,
        force = 1,
        force_pull = 1,
        max.time = 0.5,
        max.iter = 10000,
        max.overlaps = 10,
        direction = "both",
        seed = NA,
        ...
){
    if ( !(is.numeric(vec.size) && vec.size > 0) ) {
        stop("'vec.size' must be a positive number.", call. = FALSE)
    }
    if ( !(is.numeric(arrow.size) && arrow.size > 0) ) {
        stop("'arrow.size' must be a positive number.", call. = FALSE)
    }
    if ( !(is.numeric(label.size) && label.size > 0) ) {
        stop("'label.size' must be a positive number.", call. = FALSE)
    }
    if ( !.is_non_empty_string(vec.color) ) {
        stop("'vec.color' must be a non-empty string specifying a colour",
             call. = FALSE)
    }
    if ( !.is_non_empty_string(label.color) ) {
        stop("'label.color' must be a non-empty string specifying a colour",
             call. = FALSE)
    }
    if( !.is_a_bool(vec.text) ){
        stop("'vec.text' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(repel.labels) ){
        stop("'repel.labels' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Names of the ordination axes
    xvar <- colnames(vector_data)[1]
    yvar <- colnames(vector_data)[2]

    # Draw vectors from the origin to the feature coordinates
    p <- p +
        geom_segment(
            data = vector_data,
            aes(
                x = 0,
                y = 0,
                xend = .data[[xvar]],
                yend = .data[[yvar]],
                linetype = .data[["Layer"]]
            ),
            arrow = arrow(length = unit(arrow.size, "cm")),
            linewidth = vec.size,
            colour = vec.color,
            inherit.aes = FALSE
        )

    # Select the appropriate text/label geometry
    FUN <- if (repel.labels) {
        if (text.labels) geom_text_repel else geom_label_repel
    } else {
        if (text.labels) geom_text else geom_label
    }

    # Common arguments shared by all text geometries
    label_args <- list(
        data = vector_data,
        mapping = aes(
            x = .data[[xvar]],
            y = .data[[yvar]],
            label = .data[["vector_label"]]
        ),
        colour = label.color,
        size = label.size,
        inherit.aes = FALSE
    )

    # Add geometry-specific arguments
    if (repel.labels) {
        label_args <- c(label_args, list(
            min.segment.length = min.segment.length,
            box.padding = box.padding,
            point.padding = point.padding,
            force = force,
            force_pull = force_pull,
            max.time = max.time,
            max.iter = max.iter,
            max.overlaps = max.overlaps,
            direction = direction,
            seed = seed
        ))
    }

    # Add feature labels
    p <- p + do.call(FUN, label_args)

    return(p)
}

################################################################################
################## Feature covariance, variance and network ####################
#
# Each function has a method for MultiAssayExperiment (result stored in
# metadata) and SingleCellExperiment (result stored in reducedDim). The
# methods only fetch the Joint-RPCA result; the work is done by the
# .plot_jrpca_*_main() helpers below.

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCACovariance",
    signature = c(x = "MultiAssayExperiment"),
    function(x, dimred, ntop = 20, correlation = FALSE, ...){
        .check_metadata_present(dimred, x)
        res <- metadata(x)[[dimred]]
        .plot_jrpca_covariance_main(res, ntop, correlation)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCACovariance",
    signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ntop = 20, correlation = FALSE, ...){
        .check_dimred_present(dimred, x)
        res <- reducedDim(x, dimred)
        .plot_jrpca_covariance_main(res, ntop, correlation)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCAExplainedVariance",
    signature = c(x = "MultiAssayExperiment"),
    function(x, dimred, by = c("layer", "feature"), ntop = 20, ...){
        .check_metadata_present(dimred, x)
        res <- metadata(x)[[dimred]]
        .plot_jrpca_explained_variance_main(res, by, ntop)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCAExplainedVariance",
    signature = c(x = "SingleCellExperiment"),
    function(x, dimred, by = c("layer", "feature"), ntop = 20, ...){
        .check_dimred_present(dimred, x)
        res <- reducedDim(x, dimred)
        .plot_jrpca_explained_variance_main(res, by, ntop)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCANetwork",
    signature = c(x = "MultiAssayExperiment"),
    function(x, dimred, ntop = 10, correlation = TRUE, threshold = 0.9, ...){
        .check_metadata_present(dimred, x)
        res <- metadata(x)[[dimred]]
        .plot_jrpca_network_main(res, ntop, correlation, threshold, ...)
    }
)

#' @rdname plotJointRPCA
#' @export
setMethod("plotJointRPCANetwork",
    signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ntop = 10, correlation = TRUE, threshold = 0.9, ...){
        .check_dimred_present(dimred, x)
        res <- reducedDim(x, dimred)
        .plot_jrpca_network_main(res, ntop, correlation, threshold, ...)
    }
)

################################ HELP FUNCTIONS ################################
#
# The Joint-RPCA feature plots share the same structure:
#   1. .plot_jrpca_*_main() checks the input,
#   2. .get_jrpca_*_data() prepares the data for plotting, and
#   3. .plot_jrpca_*() creates the plot.
#
# All of them are based on the joint feature loadings V
# (attr(res, "rotation"); features x components) and the singular values s
# (attr(res, "varExplained"); one per component). The feature covariance is
#   cov = V diag(s)^2 V^T,
# which is how gemelli calculates it (feature_covariance_table()). The full
# matrix has one row and column per feature, so it is never calculated for all
# features. Instead, only the selected features are included.

############################ Check input and plot ##############################

# Covariance heatmap: check input, prepare long data and plot
.plot_jrpca_covariance_main <- function(res, ntop, correlation){
    .check_jrpca_input(res, ntop, correlation)
    df <- .get_jrpca_covariance_data(res, ntop, correlation)
    .plot_jrpca_covariance(df, correlation)
}

# Explained variance heatmap: check input, prepare long data and plot
.plot_jrpca_explained_variance_main <- function(res, by, ntop){
    by <- match.arg(by, c("layer", "feature"))
    .check_jrpca_input(res, ntop)
    df <- .get_jrpca_explained_variance_data(res, by, ntop)
    .plot_jrpca_explained_variance(df, by)
}

# Network: check input, prepare graph and plot
.plot_jrpca_network_main <- function(res, ntop, correlation, threshold, ...){
    .check_jrpca_input(res, ntop, correlation)
    if( !.is_a_numeric(threshold) ){
        stop("'threshold' must be a single number.", call. = FALSE)
    }
    data <- .get_jrpca_network_data(res, ntop, correlation, threshold)
    .plot_jrpca_network(data, correlation, ...)
}

# Checks shared by all Joint-RPCA feature plots
.check_jrpca_input <- function(res, ntop, correlation = FALSE){
    # The result must come from getJointRPCA() / addJointRPCA()
    if( !inherits(res, "JointRPCA") ){
        stop("The plotted results must be a class 'JointRPCA'.",
            call. = FALSE)
    }
    # ntop is the number of features per layer, NULL means all features
    if( !(is.null(ntop) || (.is_an_integer(ntop) && ntop > 0)) ){
        stop("'ntop' must be a positive integer or NULL.", call. = FALSE)
    }
    if( !.is_a_bool(correlation) ){
        stop("'correlation' must be TRUE or FALSE.", call. = FALSE)
    }
}

############################### Data preparation ###############################

# Get the data layer (e.g., "microbiota" or "metabolites") of each feature.
# Joint-RPCA stacks the features of all layers into one loading matrix, in the
# same order as the layers. attr(res, "n_features") gives the number of
# features in each layer, so the layer of each row can be recovered by
# repeating the layer names.
#' @importFrom stats setNames
.get_jrpca_layer <- function(res){
    n <- attr(res, "n_features")
    setNames(rep(names(n), times = n), rownames(attr(res, "rotation")))
}

# Variance of each feature on each component: s_j^2 V_ij^2. The result is a
# features x components matrix. Its row sums are the feature variances (the
# diagonal of the covariance matrix), and its column sums are s_j^2 because the
# loadings are orthonormal.
.get_jrpca_component_variance <- function(res){
    s <- attr(res, "varExplained")
    sweep(attr(res, "rotation")^2, 2L, s^2, "*")
}

# Select the 'ntop' features with the largest variance from each layer.
# Selecting per layer ensures that small layers (e.g., few metabolites) are not
# crowded out by large ones (e.g., many taxa). ntop = NULL keeps all features.
.select_jrpca_features <- function(res, ntop){
    # Feature variance = diagonal of the covariance matrix
    score <- rowSums(.get_jrpca_component_variance(res))
    if( is.null(ntop) ){
        return(names(score))
    }
    # Split the variances by layer and keep the largest ones within each layer
    keep <- lapply(split(score, .get_jrpca_layer(res)), function(x)
        names(sort(x, decreasing = TRUE))[seq_len(min(ntop, length(x)))])
    unlist(keep, use.names = FALSE)
}

# Covariance (or correlation) matrix of the selected features. Each element
# depends only on its own two features, so calculating it for a subset gives
# the same values as calculating the full matrix and subsetting it.
#' @importFrom stats cov2cor
.get_jrpca_matrix <- function(res, ntop, correlation){
    keep <- .select_jrpca_features(res, ntop)
    V <- attr(res, "rotation")[keep, , drop = FALSE]
    s <- attr(res, "varExplained")
    # cov = V diag(s)^2 V^T
    mat <- V %*% diag(s^2, nrow = length(s)) %*% t(V)
    # Correlation scales covariance by the feature standard deviations,
    # cor_ij = cov_ij / sqrt(cov_ii cov_jj), as in gemelli
    # feature_correlation_table().
    if( correlation ){
        mat <- cov2cor(mat)
    }
    return(mat)
}

# Convert a matrix to the long format used by .plot_heatmap(): one row per
# cell with columns FeatureID (y axis), SampleID (x axis) and value. Here
# "SampleID" holds whatever is on the columns of the matrix (a second feature
# or a component). Rows are ordered by hierarchical clustering so that similar
# features are next to each other. Layer of each row is added for facetting.
#' @importFrom stats hclust dist
.jrpca_to_long <- function(mat, layer){
    df <- as.data.frame(as.table(mat), responseName = "value")
    colnames(df)[1:2] <- c("FeatureID", "SampleID")
    # Clustering needs at least two rows
    ord <- if( nrow(mat) > 1L ) rownames(mat)[hclust(dist(mat))$order]
        else rownames(mat)
    df[["FeatureID"]] <- factor(df[["FeatureID"]], levels = ord)
    df[["Layer"]] <- layer[as.character(df[["FeatureID"]])]
    return(df)
}

# Data for the covariance heatmap. Both axes are features: FeatureID is the
# row feature and SampleID the column feature. Columns use the same clustering
# order as rows, so that the diagonal is preserved. The layers of both
# features are added so that the heatmap can be split into within- and
# between-layer blocks.
.get_jrpca_covariance_data <- function(res, ntop, correlation){
    layer <- .get_jrpca_layer(res)
    df <- .jrpca_to_long(.get_jrpca_matrix(res, ntop, correlation), layer)
    df[["SampleID"]] <- factor(
        df[["SampleID"]], levels = levels(df[["FeatureID"]]))
    df[["ColLayer"]] <- layer[as.character(df[["SampleID"]])]
    return(df)
}

# Data for the explained variance heatmap. Rows are layers (by = "layer") or
# features (by = "feature"), columns are components. Each row is scaled to sum
# to one, so the values show how the variance of a layer or feature is divided
# across the components.
.get_jrpca_explained_variance_data <- function(res, by, ntop){
    layer <- .get_jrpca_layer(res)
    # Variance of each feature on each component
    expl <- .get_jrpca_component_variance(res)
    if( by == "layer" ){
        # Sum feature variances within each layer. Each row is now a layer,
        # so the "layer" of a row is its own name (used for facetting).
        expl <- rowsum(expl, layer)
        layer <- setNames(rownames(expl), rownames(expl))
    } else{
        # Keep only the selected features
        expl <- expl[.select_jrpca_features(res, ntop), , drop = FALSE]
    }
    # Scale each row to proportions
    df <- .jrpca_to_long(expl / rowSums(expl), layer)
    return(df)
}

# Data for the network. Edges are feature pairs whose absolute covariance (or
# correlation) is at least 'threshold'. Only the upper triangle is used, so
# each pair is included once and self-loops (the diagonal) are excluded.
# plotRowGraph() takes node information from rowData, so the layers are
# stored in a small TreeSE whose rows are the nodes.
#' @importFrom tidygraph tbl_graph
#' @importFrom S4Vectors DataFrame
.get_jrpca_network_data <- function(res, ntop, correlation, threshold){
    mat <- .get_jrpca_matrix(res, ntop, correlation)
    # Row and column indices of the feature pairs that pass the threshold
    idx <- which(upper.tri(mat) & abs(mat) >= threshold, arr.ind = TRUE)
    edges <- data.frame(from = idx[, 1], to = idx[, 2], value = mat[idx])
    # Edge attributes for plotting: the sign of the relationship (edge colour)
    # and its strength (edge width)
    edges[["Sign"]] <- ifelse(edges[["value"]] > 0, "positive", "negative")
    edges[["Strength"]] <- abs(edges[["value"]])
    # Nodes are all selected features, also those without edges
    nodes <- data.frame(name = rownames(mat))
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
    se <- TreeSummarizedExperiment(rowData = DataFrame(
        Layer = .get_jrpca_layer(res)[rownames(mat)],
        row.names = rownames(mat)))
    return(list(graph = graph, se = se))
}

################################ Visualization #################################
#
# The plots reuse .plot_heatmap() and plotRowGraph() with their default colour
# scales, which are colour-blind friendly: blue-white-red for signed values,
# viridis for proportions, and blue-orange for the sign of network edges. The
# legend title comes from the name of the plotted column, so the "value"
# column is renamed before plotting.

# Covariance heatmap
.plot_jrpca_covariance <- function(df, correlation){
    name <- if( correlation ) "Correlation" else "Covariance"
    colnames(df)[colnames(df) == "value"] <- name
    # scale = TRUE gives a blue-white-red scale where zero is white
    p <- .plot_heatmap(df, fill = name, scale = TRUE)
    # Split into layer blocks. space = "free" sizes each panel by its number
    # of features, so that all cells are equal size.
    p <- p +
        facet_grid(Layer ~ ColLayer, scales = "free", space = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}

# Explained variance heatmap
.plot_jrpca_explained_variance <- function(df, by){
    name <- "Explained variance"
    colnames(df)[colnames(df) == "value"] <- name
    # scale = FALSE gives a viridis scale
    p <- .plot_heatmap(df, fill = name, scale = FALSE)
    if( by == "layer" ){
        # Few cells, so the percentages can be shown as text
        p <- p + geom_text(
            aes(label = sprintf("%.0f%%", 100 * .data[[name]])),
            colour = "white")
    } else{
        # Split features by layer
        p <- p + facet_grid(Layer ~ ., scales = "free", space = "free")
    }
    return(p)
}

# Network. Nodes are shaped by layer. Edges are coloured by the sign of the
# relationship and their width shows its strength. Nodes are not coloured, so
# that node and edge colours do not mix. Additional arguments, e.g.,
# show.label, are passed to plotRowGraph().
.plot_jrpca_network <- function(
        data, correlation, layout = "fr", edge.type = "link", ...){
    p <- plotRowGraph(
        data[["graph"]], data[["se"]], shape.by = "Layer",
        edge.colour.by = "Sign", edge.width.by = "Strength",
        layout = layout, edge.type = edge.type, ...)
    return(p)
}
