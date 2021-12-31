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
#' @param rank,... additional arguments
#' \itemize{
#'   \item{If \code{!is.null(rank)} matching arguments are passed on to
#'     \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'     \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'   }
#'   \item{additional arguments for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}}
#' }
#'   
#' @param abund_values a \code{character} value defining which assay data to
#'   use. (default: \code{abund_values = "relabundance"})
#'   
#' @param as_relative logical scalar: Should the detection threshold be applied
#'   on compositional (relative) abundances? Passed onto
#'   \code{\link[mia:getPrevalence]{getPrevalence}}. (default: \code{TRUE})
#'   
#' @param colour_by Specification of a feature to colour points by, see the 
#'   \code{by} argument in 
#'   \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for 
#'   possible values. Only used with \code{layout = "point"}.
#' @param shape_by Specification of a feature to shape points by, see the 
#'   \code{by} argument in 
#'   \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for 
#'   possible values. Only used with \code{layout = "point"}.
#' @param size_by Specification of a feature to size points by, see the 
#'   \code{by} argument in 
#'   \code{\link[scater:retrieveFeatureInfo]{?retrieveFeatureInfo}} for 
#'   possible values. Only used with \code{layout = "point"}.
#'   
#' @param facet_by Taxonomic rank to facet the plot by. 
#' Value must be of \code{taxonomyRanks(x)}
#' Argument can only be used in function plotPrevalentAbundance. 
#' 
#' @param label a \code{logical}, \code{character} or \code{integer} vector
#'   for selecting labels from the rownames of \code{x}. If \code{rank} is not 
#'   \code{NULL} the rownames might change. (default: \code{label = NULL})
#'
#' @param detections Detection thresholds for absence/presence. Either an
#'   absolutes value compared directly to the values of \code{x} or a relative
#'   value between 0 and 1, if \code{as_relative = TRUE}.
#'   
#' @param prevalences Prevalence thresholds (in 0 to 1). The
#'   required prevalence is strictly greater by default. To include the
#'   limit, set \code{include_lowest} to \code{TRUE}.
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
#' @details 
#' Agglomeration on different taxonomic levels is available through the 
#' \code{rank} argument. 
#' 
#' To exclude certain taxa, preprocess \code{x} to your liking, for example 
#' with subsetting via \code{getPrevalentTaxa} or 
#' \code{agglomerateByPrevalence}.
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
#'                    
#' # point layout for plotTaxaPrevalence can be used to visualize by additional
#' # information
#' plotPrevalentAbundance(GlobalPatterns, rank = "Family",
#'                        colour_by = "Phylum") +
#'     scale_x_log10()
#' 
#' # When using function plotPrevalentAbundace, it is possible to create facets
#' # with 'facet_by'.
#' plotPrevalentAbundance(GlobalPatterns, rank = "Family",
#'                        colour_by = "Phylum", facet_by = "Kingdom") +
#'     scale_x_log10()
NULL

################################################################################
# plotPrevalence

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
        .check_assay_present(abund_values, x)
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
                                               BPPARAM)
        plot_data$colour_by <- plot_data$colour_by * 100
        .prevalence_plotter(plot_data, 
                            layout = "line",
                            xlab = ifelse(as_relative,"Abundance [%]","Detection"),
                            ylab = "N",
                            colour_by = "Prevalence [%]",
                            size_by = NULL,
                            shape_by = NULL,
                            ...)
    }
)

#' @importFrom DelayedArray rowMeans
#' @importFrom BiocGenerics ncol
.get_prevalence_value <- function(d, mat){
    rowSums(mat > d) / ncol(mat)
}

.get_prevalence_count <- function(d, p, mat){
    sum(.get_prevalence_value(d, mat) >= p)
}

#' @importFrom BiocParallel bpmapply bpisup bpstart bpstop SerialParam
#' @importFrom SummarizedExperiment assay
.get_prevalence_plot_data <- function(x, abund_values, detections, prevalences,
                                      as_relative = TRUE, 
                                      BPPARAM = BiocParallel::SerialParam()){
    mat <- assay(x, abund_values, withDimnames = TRUE)
    if(as_relative){
        mat <- mia:::.calc_rel_abund(mat)
    }
    ans <- expand.grid(detection = detections, prevalence = prevalences)
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    ans$n <- unlist(bpmapply(.get_prevalence_count,
                             ans$detection,
                             ans$prevalence,
                             MoreArgs = list(mat = mat),
                             BPPARAM = BPPARAM,
                             SIMPLIFY = FALSE))
    colnames(ans) <- c("X","colour_by","Y")
    ans
}

################################################################################
# plotPrevalentAbundance

#' @rdname plotPrevalence
#' @export
setGeneric("plotPrevalentAbundance", signature = c("x"),
           function(x, ...) standardGeneric("plotPrevalentAbundance"))

#' @rdname plotPrevalence
#' @export
setMethod("plotPrevalentAbundance", signature = c(x = "SummarizedExperiment"),
    function(x,
             rank = taxonomyRanks(x)[1L],
             abund_values = "counts",
             as_relative = TRUE,
             colour_by = NULL,
             size_by = NULL,
             shape_by = NULL,
             label = NULL,
             facet_by = NULL,
             ...){
        # input check
        .check_assay_present(abund_values, x)
        if(!.is_a_bool(as_relative)){
            stop("'as_relative' must be TRUE or FALSE.", call. = FALSE)
        }
        # Check facet_by. It is FALSE by default, but user can specify it, but 
        # the value must be in taxonomyRanks.
        if(!(is.null(facet_by) || facet_by %in% taxonomyRanks(x))){
            stop("'facet_by' must be in taxonomyRanks.",  call. = FALSE)
        }
        
        x <- mia:::.agg_for_prevalence(x, rank, na.rm = TRUE, relabel = TRUE,
                                       ...)
        label <- .norm_label(label, x)
        #
        plot_data <- .get_prevalence_plot_point_data(x, abund_values, 
                                                     as_relative = as_relative,
                                                     label = label)
        vis_out <- .incorporate_prevalence_vis(plot_data,
                                               se = x,
                                               colour_by = colour_by,
                                               size_by = size_by,
                                               shape_by = shape_by,
                                               label = label,
                                               facet_by = facet_by)
        plot_data <- vis_out$df
        colour_by <- vis_out$colour_by
        size_by <- vis_out$size_by
        shape_by <- vis_out$shape_by
        facet_by <- vis_out$facet_by
        xlab <- paste0(ifelse(as_relative, "Rel. ", ""),"Abundance")
        ylab <- paste0("Prevalence(", ifelse(is.null(rank), "Features", rank),
                       ") [%]")
        plot <- .prevalence_plotter(plot_data, 
                            layout = "point",
                            xlab = xlab,
                            ylab = ylab,
                            colour_by = colour_by,
                            size_by = size_by,
                            shape_by = shape_by,
                            ...)

        # If facet_by is not NULL, user has specified it. Adds the facets to the plot.
        if(!is.null(facet_by)){
            plot <- plot + 
                # Create facets
                facet_wrap(vars(!!sym("facet_by")))
        }
        
        return(plot)
    }
)

#' @importFrom DelayedArray rowMeans
#' @importFrom SummarizedExperiment assay
.get_prevalence_plot_point_data <- function(x, abund_values, as_relative = TRUE,
                                            label = NULL){
    mat <- assay(x, abund_values, withDimnames = TRUE)
    if(as_relative){
        mat <- mia:::.calc_rel_abund(mat)
    }
    ans <- data.frame(X = rowMeans(mat, na.rm = TRUE),
                      Y = .get_prevalence_value(0, mat) * 100,
                      label = rownames(mat))
    if(!is.null(label)){
        ans$label[!label] <- NA
    }
    ans
}

#' @importFrom scater retrieveFeatureInfo
.incorporate_prevalence_vis <- function(plot_data,
                                        se = se,
                                        colour_by = NULL,
                                        size_by = NULL,
                                        shape_by = NULL,
                                        label = NULL,
                                        facet_by = NULL){
    variables <- c(colour_by = colour_by,
                   size_by = size_by,
                   shape_by = shape_by,
                   facet_by = facet_by)
    if(!is.null(variables)){
        for(i in seq_along(variables)){
            # get data
            feature_info <- retrieveFeatureInfo(se, variables[i],
                                                search = "rowData")
            # mirror back variable name, if a partial match was used
            var_name <- names(variables)[i]
            assign(var_name, 
                   .get_new_var_name_value(get(var_name),
                                           feature_info$name))
            plot_data[,names(feature_info$name)] <- feature_info$value
            if(!is.null(label)){
                if(is.factor(feature_info$value) || 
                   is.character(feature_info$value)){
                    plot_data[,names(feature_info$name)][!label] <- NA
                }
            }
        }
    }
    return(list(df = plot_data,
                colour_by = colour_by,
                size_by = size_by,
                shape_by = shape_by,
                facet_by = facet_by))
}

################################################################################
# plotTaxaPrevalence

#' @rdname plotPrevalence
#' @export
setGeneric("plotTaxaPrevalence", signature = c("x"),
           function(x, ...) standardGeneric("plotTaxaPrevalence"))

#' @rdname plotPrevalence
#' @export
setMethod("plotTaxaPrevalence", signature = c(x = "SummarizedExperiment"),
          function(x,
                   rank = taxonomyRanks(x)[1L],
                   abund_values = "counts",
                   detections = NULL,
                   ndetections = 20,
                   as_relative = TRUE,
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
        .check_assay_present(abund_values, x)
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
        x <- mia:::.agg_for_prevalence(x, rank, na.rm = TRUE, relabel = TRUE,
                                       ...)
        plot_data <- .get_prevalence_plot_matrix(x, abund_values, detections,
                                                 as_relative, 
                                                 min_prevalence,
                                                 BPPARAM)
        plot_data$colour_by <- plot_data$colour_by * 100
        xlab <- ifelse(as_relative,"Abundance [%]","Detection")
        ylab <- ifelse(is.null(rank), "Features", rank)
        colour_by <- "Prevalence [%]"
        .prevalence_plotter(plot_data, 
                            layout = "heatmap",
                            xlab = xlab,
                            ylab = ylab,
                            colour_by = colour_by,
                            size_by = NULL,
                            shape_by = NULL,
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
    mat <- assay(x, abund_values, withDimnames = TRUE)
    if(as_relative){
        mat <- mia:::.calc_rel_abund(mat)
    }
    
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    ans <- bplapply(detections,
                    .get_prevalence_value,
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
                                layout = c("line","point","heatmap"),
                                xlab = NULL,
                                ylab = NULL,
                                colour_by = NULL,
                                size_by = NULL,
                                shape_by = NULL,
                                flipped = FALSE,
                                add_legend = TRUE,
                                point_alpha = 1,
                                point_size = 2,
                                line_alpha = 1,
                                line_type = NULL,
                                line_size = 1){
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
    } else if(layout == "point"){
        point_args <- .get_point_args(colour_by = colour_by, shape_by = shape_by,
                                      size_by = size_by,
                                      alpha = point_alpha,
                                      size = point_size)
        plot_out <- plot_out +
            do.call(geom_point, point_args$args)
        # resolve the colours
        plot_out <- .resolve_plot_colours(plot_out,
                                          plot_data$colour_by,
                                          colour_by,
                                          fill = TRUE,
                                          na.translate = FALSE)
        # add additional guides
        plot_out <- .add_extra_guide(plot_out, shape_by, size_by)
    } else if(layout == "heatmap"){
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
    } else {
        stop("incompatible layout")
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
