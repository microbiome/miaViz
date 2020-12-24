#' Plotting abundance data
#'
#' \code{plotAbundance} plots the abundance on a selected taxonomic rank.
#' Since this probably makes sense only for relative abundance data, the
#' assay used by default is expected to be in the slot \sQuote{relabundance}.
#'
#' Subsetting to rows of interested and ordering of those is expected to be done
#' outside of this functions, e.g. \code{x[1:2,]}. This will plot data of all
#' features present.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'
#' @param rank a single \code{character} value defining the taxonomic rank to
#'   use. Must be a value of \code{taxonomyRanks(x)}.
#'
#' @param abund_values a \code{character} value defining which assay data to
#'   use. (default: \code{abund_values = "relabundance"})
#'
#' @param layout Either \dQuote{bar} or \dQuote{point}. 
#' 
#' @param one_facet Should the plot be returned in on facet or split into 
#'   different facet, one facet per different value detect in \code{rank}
#' 
#' @param ncol,scales if \code{one_facet = FALSE}, \code{ncol} defines many 
#'   columns should be for plotting the different facets and \code{scales} is
#'   used to define the behavior of the scales of each facet. Both values are 
#'   passed onto \code{\link[ggplot2:facet_wrap]{facet_wrap}}.
#' 
#' @param ... additional parameters for plotting
#'
#' @return a \code{\link[ggplot2:ggplot]{ggplot}} object
#'
#' @name plotAbundance
#'
#' @examples
#' data(GlobalPatterns)
#' se <- GlobalPatterns
#' 
#' #
#' plotAbundance(se, abund_values="counts")
#' #
#' plotAbundance(se, abund_values="counts", rank = "Phylum", add_legend = FALSE)
#' 
#' # If rank is set to NULL plotAbundance behaves like plotExpression
#' plotAbundance(se, abund_values="counts", rank = NULL,
#'               features = head(rownames(se)))
#' 
NULL

#' @rdname plotAbundance
setGeneric("plotAbundance", signature = c("x"),
           function(x, ...)
               standardGeneric("plotAbundance"))

.check_abund_plot_args <- function(one_facet = TRUE,
                                   ncol = 2){
    
   if(!.is_a_bool(one_facet)){
        stop("'one_facet' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!is.numeric(ncol) || as.integer(ncol) != ncol || ncol < 1){
        stop("'ncol' must be an integer value above or equal to 1.",
             call. = FALSE)
    }
}


#' @rdname plotAbundance
#' @importFrom scater plotExpression
#' @importFrom ggplot2 facet_wrap
#' @export
setMethod("plotAbundance", signature = c("SummarizedExperiment"),
    function(x,
             rank = taxonomyRanks(x)[1],
             abund_values = "relabundance",
             layout = c("bar","point"),
             one_facet = TRUE,
             ncol = 2,
             scales = "fixed",
             ...){
        # input checks
        if(nrow(x) == 0L){
            stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        # if rank is set to NULL, default to plotExpression 
        if(is.null(rank)){
            plot <- plotExpression(x, exprs_values = abund_values, one_facet = one_facet,
                                   ncol = ncol, scales = scales, ...)
            ylab <- gsub("Expression","Abundance",plot$labels$y)
            plot <- plot +
                ylab(ylab)
            return(plot)
        }
        # input checks
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be an non empty single character value.",
                 call. = FALSE)
        }
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        layout <- match.arg(layout, c("bar","point"))
        .check_abund_plot_args(one_facet = one_facet,
                               ncol = ncol)
        #
        abund_data <- .get_abundance_data(x, rank, abund_values)
        plot_out <- .abund_plotter(abund_data,
                                   xlab = "Samples",
                                   ylab = "Rel. Abundance",
                                   colour_by = rank,
                                   layout = layout,
                                   ...)
        if (!one_facet) {
            plot_out <- plot_out + 
                facet_wrap(~colour_by, ncol = ncol, scales = scales)
        }
        plot_out
    }
)


MELT_NAME <- "Sample"
MELT_VALUES <- "Value"

#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr mutate group_by summarize rename
#' @importFrom tidyr pivot_longer nest unnest
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
.get_abundance_data <- function(x, rank, abund_values){
    data <- assay(x,abund_values)
    data <- sweep(data, 2L, colSums(data),"/")
    merge_FUN <- function(data){
        data %>%
            group_by(!!sym(MELT_NAME)) %>%
            summarize(!!sym(MELT_VALUES) := sum(!!sym(MELT_VALUES)))
    }
    data <- data %>%
        as.data.frame() %>%
        mutate(rank = factor(rowData(x)[,rank], unique(rowData(x)[,rank]))) %>%
        pivot_longer(cols = !rank,
                     names_to = MELT_NAME,
                     values_to = MELT_VALUES) %>%
        mutate(!!sym(MELT_NAME) := factor(!!sym(MELT_NAME), unique(!!sym(MELT_NAME)))) %>%
        nest(data = !rank) %>%
        mutate(data = map(data, merge_FUN)) %>%
        unnest(cols = data)
    # rename columns
    data <- data %>%
        dplyr::rename(colour_by = "rank",
                      X = MELT_NAME,
                      Y = MELT_VALUES)
    data
}


################################################################################
# abundance plotter


#' @importFrom ggplot2 ggplot theme_classic geom_point geom_bar coord_flip
#'   scale_y_continuous
.abund_plotter <- function(object,
                           xlab = NULL,
                           ylab = NULL,
                           colour_by = NULL,
                           layout = "bar",
                           flipped = FALSE,
                           add_legend = TRUE,
                           add_x_text = FALSE,
                           bar_alpha = 0.65,
                           point_alpha = 1,
                           point_size = 2){
    # start plotting
    plot_out <- ggplot(object, aes_string(x="X", y="Y")) +
        xlab(xlab) +
        ylab(ylab)
    # either bar or point plot
    if(layout == "bar"){
        abund_out <- .get_bar_args(colour_by,
                                   colour_by,
                                   alpha = bar_alpha)
        plot_out <- plot_out +
            do.call(geom_bar, c(abund_out$args, list(stat="identity"))) +
            scale_y_continuous(expand = c(0,0))
        
    } else {
        abund_out <- .get_point_args(colour_by,
                                     shape_by = NULL,
                                     size_by = NULL,
                                     alpha = point_alpha,
                                     size = point_size)
        plot_out <- plot_out +
            do.call(geom_point, abund_out$args)
    }
    # adjust point colours
    if(!is.null(colour_by)){
        # resolve the colour for the line colours
        plot_out <- .resolve_plot_colours(plot_out,
                                          object$colour_by,
                                          colour_by,
                                          fill = FALSE)
        # resolve the color for fill
        plot_out <- .resolve_plot_colours(plot_out,
                                          object$colour_by,
                                          colour_by,
                                          fill = TRUE)
    }
    plot_out <- plot_out +
        theme_classic()
    # add legend
    if(!add_legend){
        plot_out <- plot_out +
            theme(legend.position = "none")
    }
    
    # flip
    if (flipped) {
        plot_out <- plot_out + coord_flip()
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
        }
    } else {
        if(!add_x_text){
            plot_out <- plot_out +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
        } else {
            plot_out <- plot_out +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
    }
    plot_out
}
