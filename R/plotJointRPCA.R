#' @name
#' plotJointRPCA
#'
#' @title
#' Visualize Joint-RPCA results
#'
#' @description
#' Creates a two-dimensional ordination plot from Joint-RPCA results stored in a
#' \code{SingleCellExperiment} or \code{MultiAssayExperiment}. In addition to the
#' sample ordination, feature loadings can be visualized as vectors.
#'
#' @details
#' This function is a wrapper around \code{\link[=plotOrdination]{plotOrdination}}
#' for Joint-RPCA results. Consequently, most graphical parameters are passed
#' directly to \code{plotOrdination()}, including options for colouring,
#' grouping, faceting and adding ellipses. See
#' \code{\link[=plotOrdination]{plotOrdination}} for a complete description of
#' these arguments.
#'
#' Feature loadings are plotted as vectors originating from the origin. By
#' default, only the longest loading vectors are shown for each data layer to
#' improve readability.
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
#' @param add.vectors: \code{logical(1)} or \code{character}. If
#' \code{TRUE}, feature loading vectors are added. Alternatively, a character
#' vector can be supplied to display only features whose names match the given
#' pattern(s). (Default: \code{TRUE})
#'
#' @param ...
#' Additional arguments passed to
#' \code{\link[=plotOrdination]{plotOrdination}} and to the feature-vector
#' plotting.
#' Commonly used Joint-RPCA-specific arguments include:
#' \itemize{
#' \item \code{ntop}: \code{integer(1)} or \code{NULL}. Maximum number of
#' loading vectors shown per data layer. The vectors with the largest lengths
#' are retained. Set to \code{NULL} to display all vectors. (Default: \code{10})
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
#' # Show all loading vectors
#' plotJointRPCA(
#'     mae,
#'     "JointRPCA",
#'     ntop = NULL
#' )
#'
#' # Display only selected feature vectors
#' plotJointRPCA(
#'     mae,
#'     "JointRPCA",
#'     add.vectors = c("Bacteroides", "Roseburia")
#' )
#'
#' # Use boxed labels without repelling
#' plotJointRPCA(
#'     mae,
#'     "JointRPCA",
#'     text.labels = FALSE,
#'     repel.labels = FALSE
#' )
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[=plotOrdination]{plotOrdination}}
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
    function(x, dimred, add.vectors = TRUE, ...){
        .check_dimred_present(dimred, x)
        if( !inherits(reducedDim(x, dimred), "JointRPCA") ){
            stop("The plotted results must be a class 'JointRPCA'.",
                call. = FALSE)
        }
        if( !( .is_a_bool(add.vectors) || is.character(add.vectors)) ){
            stop("'add.vectors must be TRUE or FALSE or character vector.",
                 call. = FALSE)
        }
        p <- plotOrdination(x, dimred, ...)
        vector_data <- .get_joint_rpca_vector_data(x, dimred, add.vectors, ...)
        if( !is.null(vector_data ) ){
            p <- .add_ordination_vectors(p, vector_data, ...)
        }
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# Get feature loadings that will be plotted as vectors
#' @importFrom dplyr group_by slice_max ungroup
.get_joint_rpca_vector_data <- function(
        tse, reduced_dim, add.vectors, ignore.case = FALSE, ntop = 10,
        ...){
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
