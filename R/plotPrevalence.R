#' Plot prevalence information
#' 
#' \code{plotPrevalence} and \code{plotTaxaPrevalence} visualize prevalence 
#' information.
#' 
#' Whereas \code{plotPrevalence} procudes a line plot, \code{plotTaxaPrevalence}
#' returns a heatmap. 
#' 
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'
#' @param detections Detection thresholds for absence/presence. Either an
#'   absolutes value compared directly to the values of \code{x} or a relative
#'   value between 0 and 1, if \code{as_relative = TRUE}.
#'   
#' @param prevalences Prevalence thresholds (in 0 to 1). The
#'   required prevalence is strictly greater by default. To include the
#'   limit, set \code{include_lowest} to \code{TRUE}.
#'   
#' @param abund_values a \code{character} value defining which assay data to
#'   use. (default: \code{abund_values = "relabundance"})
#'   
#' @param as_relative logical scalar: Should the detection threshold be applied
#'   on compositional (relative) abundances? Passed onto
#'   \code{\link[mia:getPrevalence]{getPrevalence}}. (default: \code{TRUE})
#' 
#' @param rank,... additional arguments
#' \itemize{
#'   \item{If \code{!is.null(rank)} matching arguments are passed on to
#'     \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'     \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#'   \item{additional arguments for plotting}
#' }
#' 
#' @param min_prevalence a single numeric value to apply as a threshold for 
#'   plotting. The threshold is applied per row and column.
#'   (default: \code{min_prevalence = 0})
#' 
#' @param ndetections If \code{detections} is \code{NULL}, a number of breaks 
#'   are calculated automatically. \code{as_relative} is then also regarded as 
#'   \code{TRUE}.
#' 
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether the UniFrac calculation should be parallelized.
#' 
#' @return 
#' A \code{ggplot2} object or \code{plotly} object, if more than one 
#' \code{prevalences} was defined.
#' 
#' @seealso 
#' \code{\link[mia:getPrevalence]{getPrevalence}},
#' \code{\link[mia:getPrevalence]{agglomerateByPrevalence}},
#' \code{\link[mia:agglomerate-methods]{agglomerateByRank}}
#' 
#' @name plotPrevalence
#' 
#' @details 
#' Agglomeration on different taxonomic levels is available through the 
#' \code{rank} argument. 
#' 
#' To exclude certain taxa, preprocess \code{x} to your likeing, for example 
#' with subsetting via \code{getPrevalentTaxa} or 
#' \code{agglomerateByPrevalence}.
#' 
#' @examples 
#' data(GlobalPatterns, package = "mia")
#' 
#' # plotting N of prevalence exceeding taxa on the Phylum level
#' plotPrevalence(GlobalPatterns, rank = "Phylum")
#' plotPrevalence(GlobalPatterns, rank = "Phylum") + scale_x_log10()
#' 
#' # plotting prevalence per taxa for different detectio thresholds as heatmap
#' plotTaxaPrevalence(GlobalPatterns, rank = "Phylum")
#' 
#' # by default a continous scale is used for different detections levels, 
#' # but this can be adjusted
#' plotTaxaPrevalence(GlobalPatterns, rank = "Phylum",
#'                    detections = c(0, 0.001, 0.01, 0.1, 0.2))
NULL

#' @rdname plotPrevalence
#' @export
setGeneric("plotPrevalence", signature = c("x"),
           function(x, ...) standardGeneric("plotPrevalence"))

#' @rdname plotPrevalence
#' @export
setMethod("plotPrevalence", signature = c(x = "SummarizedExperiment"),
          function(x,
                   detections = c(0.01, 0.1, 1, 2, 5, 10, 20)/100,
                   prevalences = seq(0.1, 1, 0.1),
                   abund_values = "counts",
                   as_relative = TRUE,
                   rank = NULL,
                   BPPARAM = BiocParallel::SerialParam(),
                   ...){
        # input check
        if(!all(.is_numeric_string(detections))){
            stop("'detections' must be numeric values.", call. = FALSE)
        }
        if(!all(.is_numeric_string(prevalences)) || any(prevalences < 0) ||
           any(prevalences > 1)){
            stop("'prevalences' must be numeric values between 0 and 1.",
                 call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        if(!.is_a_bool(as_relative)){
            stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
        }
        if(as_relative && (any(detections < 0) || any(detections > 1))){
            stop("If 'as_relative' == TRUE, detections' must be numeric ",
                 "values between 0 and 1.", call. = FALSE)
        }
        #
        x <- mia:::.agg_for_prevalence(x, rank, ...)
        plot_data <- .get_prevalence_plot_data(x, abund_values, detections,
                                               prevalences, as_relative,
                                               min_prevalence,
                                               BPPARAM)
        plot_data$colour_by <- plot_data$colour_by * 100
        .prevalence_plotter(plot_data, 
                            layout = "line",
                            xlab = paste0("Detection ",
                                          ifelse(as_relative," [%]","")),
                            ylab = "N",
                            colour_by = "Prevalence [%]",
                            ...)
    }
)

#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
#' @importFrom SummarizedExperiment assay
.get_prevalence_plot_data <- function(x, abund_values, detections, prevalences,
                                      as_relative = TRUE, 
                                      BPPARAM = BiocParallel::SerialParam()){
    mat <- assay(x,abund_values)
    if(as_relative){
        mat <- mia:::.calc_rel_abund(mat)
    }
    FUN_prevalences <- function(d, p, mat){
        sum((rowSums(mat > d) / ncol(mat)) >= p)
    }
    ans <- expand.grid(detection = detections, prevalence = prevalences)
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    ans$n <- unlist(bpmapply(FUN_prevalences,
                             ans$detection,
                             ans$prevalence,
                             MoreArgs = list(mat = mat),
                             BPPARAM = BPPARAM,
                             SIMPLIFY = FALSE))
    colnames(ans) <- c("X","colour_by","Y")
    ans
}


#' @rdname plotPrevalence
#' @export
setGeneric("plotTaxaPrevalence", signature = c("x"),
           function(x, ...) standardGeneric("plotTaxaPrevalence"))

#' @rdname plotPrevalence
#' @export
setMethod("plotTaxaPrevalence", signature = c(x = "SummarizedExperiment"),
          function(x,
                   detections = NULL,
                   ndetections = 20,
                   abund_values = "counts",
                   as_relative = TRUE,
                   rank = taxonomyRanks(x)[1L],
                   min_prevalence = 0,
                   BPPARAM = BiocParallel::SerialParam(),
                   ...){
        # input check
        if(!is.null(detections)){
            if(!all(.is_numeric_string(detections))){
                stop("'detections' must be numeric values.", call. = FALSE)
            }
        } else {
            if(!is.numeric(ndetections) || 
               ndetections != as.integer(ndetections) ||
               length(ndetections) != 1L){
                stop("'ndetections' must be a single integer value.",
                     call. = FALSE)
            }
            detections <- seq(0,1,length.out = ndetections + 1L)
            as_relative <- TRUE
        }
        .check_abund_values(abund_values, x)
        if(!.is_a_bool(as_relative)){
            stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
        }
        if(as_relative && (any(detections < 0) || any(detections > 1))){
            stop("If 'as_relative' == TRUE, detections' must be numeric ",
                 "values between 0 and 1.", call. = FALSE)
        }
        if(length(min_prevalence) != 1 || !.is_numeric_string(min_prevalence)){
            stop("'min_prevalence' must be single numeric values.",
                 call. = FALSE)
        }
        #
        x <- mia:::.agg_for_prevalence(x, rank, na.rm = TRUE, ...)
        rownames(x) <- getTaxonomyLabels(x, make_unique = TRUE) 
        plot_data <- .get_prevalence_plot_matrix(x, abund_values, detections,
                                                 as_relative, 
                                                 min_prevalence,
                                                 BPPARAM)
        plot_data$colour_by <- plot_data$colour_by * 100
        .prevalence_plotter(plot_data, 
                            layout = "heatmap",
                            xlab = paste0("Detection ",
                                          ifelse(as_relative," [%]","")),
                            ylab = ifelse(is.null(rank), "Features", rank),
                            colour_by = "Prevalence [%]",
                            ...)
    }
)

.is_continuous <- function(i){
    i <- sort(unique(i))
    z <- round(c(i[-1L],max(i)) - i,5L)
    length(unique(z[-length(z)])) == 1L
}

#' @importFrom BiocParallel bplapply bpisup bpstart bpstop SerialParam
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#' @importFrom DelayedArray rowSums
.get_prevalence_plot_matrix <- function(x, abund_values, detections, 
                                        as_relative = TRUE, 
                                        min_prevalence,
                                        BPPARAM = BiocParallel::SerialParam()){
    mat <- assay(x,abund_values)
    if(as_relative){
        mat <- mia:::.calc_rel_abund(mat)
    }
    FUN_prevalences <- function(d, mat){
        rowSums(mat > d) / ncol(mat)
    }
    
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    ans <- bplapply(detections,
                    FUN_prevalences,
                    mat = mat,
                    BPPARAM = BPPARAM)
    ans <- data.frame(ans)
    colnames(ans) <- detections
    f <- ans >= min_prevalence
    ans <- ans[rowSums(f) != 0,colSums(f) != 0,drop=FALSE]
    if(any(dim(ans) == 0)){
        stop("No data left after apply threshold 'min_prevalence'.",
             call. = FALSE)
    }
    lvls <- rownames(ans)[order(rowSums(ans))]
    ans$ID <- rownames(mat)[rowSums(f) != 0]
    ans <- ans %>%
        pivot_longer(!.data$ID, 
                     names_to = "detection",
                     values_to = "prevalence")
    colnames(ans) <- c("Y","X","colour_by")
    ans$X <- round(as.numeric(ans$X),4) * 100
    if(!.is_continuous(ans$X)){
        ans$X <- factor(ans$X,
                        unique(as.character(sort(as.numeric(ans$X)))))
    }
    ans$Y <- factor(ans$Y,lvls)
    ans
}

#' @importFrom ggplot2 ggplot geom_raster geom_point geom_line
#'   scale_fill_distiller scale_y_discrete scale_x_discrete scale_x_continuous
#'   theme_classic
.prevalence_plotter <- function(plot_data,
                                layout = c("line","heatmap"),
                                xlab = NULL,
                                ylab = NULL,
                                colour_by = NULL,
                                flipped = FALSE,
                                add_legend = TRUE,
                                point_alpha = 1,
                                point_size = 2,
                                line_alpha = 1,
                                line_type = NULL,
                                line_size = 1,
                                ...){
    plot_out <- ggplot(plot_data, aes_string(x = "X", y = "Y")) +
        labs(x = xlab, y = ylab)
    if(layout == "line"){
        point_args <- .get_point_args(colour_by = colour_by, shape_by = NULL,
                                      size_by = NULL,
                                      alpha = point_alpha,
                                      size = point_size)
        line_args <- .get_line_args(colour_by = colour_by, linetype_by = NULL,
                                    size_by = NULL,
                                    alpha = line_alpha,
                                    linetype = line_type,
                                    size = line_size)
        point_args$args$mapping$group <- sym("colour_by")
        line_args$args$mapping$group <- sym("colour_by")
        plot_out <- plot_out +
            do.call(geom_point, point_args$args) +
            do.call(geom_line, line_args$args)
        # resolve the colours
        plot_out <- .resolve_plot_colours(plot_out,
                                          plot_data$colour_by,
                                          colour_by,
                                          fill = TRUE)
        plot_out <- .resolve_plot_colours(plot_out,
                                          plot_data$colour_by,
                                          colour_by,
                                          fill = FALSE)
    } else {
        raster_args <- .get_bar_args(colour_by = colour_by, alpha = 1,
                                     add_border = FALSE)
        plot_out <- plot_out +
            do.call(geom_raster, raster_args$args) +
            scale_fill_distiller(palette = "RdYlBu", name = colour_by) +
            scale_y_discrete(expand = c(0,0))
        if(is.factor(plot_data$X)){
            plot_out <- plot_out + 
                scale_x_discrete(expand = c(0,0))
        } else {
            plot_out <- plot_out + 
                scale_x_continuous(expand = c(0,0), n.breaks = 7L)
        }
    }
    # theme
    plot_out <- plot_out +
        theme_classic()
    # add legend
    plot_out <- .add_legend(plot_out, add_legend)
    # flip
    plot_out <- .flip_plot(plot_out, flipped, add_x_text = TRUE,
                           angle_x_text = FALSE)
    plot_out
}
