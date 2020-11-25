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
#' @param across,per a \code{character} value selecting the values to use
#'   from \code{colData(x)}. For \code{across} a continous numeric value or a
#'   factor is expected and for \code{per} a factor or character value. By
#'   default the first occurance of these types is used for plotting. (default:
#'   \code{across = sample_val = NULL})
#'
#' @param abund_values a \code{character} value defining which assay data to
#'   use. (default: \code{abund_values = "relabundance"})
#'
#' @return a \code{\link[ggplot2:ggplot]{ggplot}} object
#'
#' @name plotAbundance
#' @aliases plotAbundanceAcross plotAbundancePer
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' plotRelAbundance(GlobalPatterns,
#'                  rank = "Phylum")
#' plotRelAbundance(GlobalPatterns,
#'                  rank = "Phylum",
#'                  subset_by = c("Crenarchaeota","Firmicutes","Bacteroidetes"),
#'                  order_by = c("Bacteroidetes","Firmicutes"))
#' }
NULL


.get_first_numeric_or_factor <- function(x){
    f <- vapply(colData(x),is.numeric,logical(1)) |
        vapply(colData(x),is.factor,logical(1))
    if(!any(f)){
        stop("No numeric of factor values found in colData(x).", call. = FALSE)
    }
    colnames(colData(x))[f][1L]
}

.get_first_factor_or_character <- function(x){
    f <- vapply(colData(x),is.character,logical(1)) |
        vapply(colData(x),is.factor,logical(1))
    if(!any(f)){
        stop("No character of factor values found in colData(x).",
             call. = FALSE)
    }
    colnames(colData(x))[f][1L]
}

.norm_across <- function(across, x){
    if(!is.null(across)){
        if(!is.character(across) || length(across) == 0L){
            stop("'across' must be a single character value.", call. = FALSE)
        }
        if(!(across %in% colnames(colData(x)))){
            stop("'across' must define a column name of colData(x)",
                 call. = FALSE)
        }
        values <- colData(x)[,across]
        if(!is.numeric(values) && !is.factor(values)){
            stop("values defined by 'across' must be numeric or a factor",
                 call. = FALSE)
        }
    } else {
        across <- .get_first_numeric_or_factor(x)
    }
    across
}

.norm_per <- function(per, x){
    if(!is.null(per)){
        if(!is.character(per) || length(per) == 0L){
            stop("'per' must be a single character value.", call. = FALSE)
        }
        if(!(per %in% colnames(colData(x)))){
            stop("'per' must define a column name of colData(x)",
                 call. = FALSE)
        }
        values <- colData(x)[,per]
        if(!is.character(values) && !is.factor(values)){
            stop("values defined by 'across' must be character or a factor",
                 call. = FALSE)
        }
    } else {
        per <- .get_first_factor_or_character(x)
    }
    per
}

#' @rdname plotAbundance
setGeneric("plotAbundance", signature = c("x"),
           function(x, ...)
               standardGeneric("plotAbundance"))

#' @rdname plotAbundance
setGeneric("plotAbundanceAcross", signature = c("x"),
           function(x, rank = taxonomyRanks(x)[1], across = NULL,
                    abund_values = "relabundance",
                    ...)
               standardGeneric("plotAbundanceAcross"))

#' @rdname plotAbundance
setGeneric("plotAbundancePer", signature = c("x"),
           function(x, rank = taxonomyRanks(x)[1], per = NULL,
                    abund_values = "relabundance",
                    ...)
               standardGeneric("plotAbundancePer"))

#' @rdname plotAbundance
setGeneric("plotAbundanceOrdered", signature = c("x"),
           function(x, ...)
               standardGeneric("plotAbundanceOrdered"))

#' @rdname plotAbundance
#' @export
setMethod("plotAbundance", signature = c("SummarizedExperiment"),
    function(x, ...){
        plot <- plotExpression(x, ...)
        ylab <- gsub("Expression","Abundance",plot$labels$y)
        plot <- plot +
          ylab(ylab)
        plot
    }
)


#' @rdname plotAbundance
#' @export
setMethod("plotAbundanceAcross", signature = c("SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], across = NULL,
             abund_values = "relabundance", scale = FALSE,
             type = c("bar","point"), legend = TRUE, ...){
        # input checks
        if(!.is_non_empty_string(rank)){
          stop("'rank' must be an non empty single character value.",
               call. = FALSE)
        }
        if(nrow(x) == 0L){
          stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        across <- .norm_across(across, x)
        if(!.is_a_bool(scale)){
          stop("'scale' must be TRUE or FALSE.", call. = FALSE)
        }
        type <- match.arg(type, c("bar","point"))
        if(!.is_a_bool(legend)){
          stop("'legend' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        plot_args <- list(theme = theme_bw)
        plot_data <- .get_abundance_plot_data(x, rank, across, abund_values,
                                            scale)
        plot_data <- .order_plot_data(plot_data, across)
        #
        plot <- .get_base_plot(plot_data, plot_args)
        if(type == "point" ||
         (!scale && length(unique(plot_data[[MELT_NAME]])) == 1L)){
          plot <- .get_rel_abundance_point_plot_across(plot, rank, across)
        } else {
          plot <- .get_rel_abundance_plot_across(plot, rank, across)
        }
        plot <- .add_abundance_plot_labels(plot, rank, across, scale,
                                         abund_values)
        plot <- .add_abundance_legend(plot, legend)
        plot
    }
)



#' @rdname plotAbundance
#' @export
setMethod("plotAbundancePer", signature = c("SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], per = NULL,
             abund_values = "relabundance", scale = FALSE, ...){
        # input checks
        if(!.is_non_empty_string(rank)){
          stop("'rank' must be an non empty single character value.",
               call. = FALSE)
        }
        if(nrow(x) == 0L){
          stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        per <- .norm_per(per, x)
        if(!.is_a_bool(scale)){
          stop("'scale' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(legend)){
          stop("'legend' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        plot_args <- list(theme = theme_bw)
        plot_data <- .get_abundance_plot_data(x, rank, per, abund_values,
                                            scale,
                                            plot_args)
        plot_data <- .order_plot_data(plot_data, per)
        #
        plot <- .get_base_plot(plot_data, plot_args)
        plot <- .get_rel_abundance_plot_per(plot, rank, per)
        plot <- .add_abundance_plot_labels(plot, rank, per, scale,
                                         abund_values)
        plot <- .add_abundance_legend(plot, legend)
        plot
    }
)

MELT_NAME <- "Sample"
MELT_VALUES <- "Value"

#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr %>% mutate left_join group_by summarize
#' @importFrom rlang sym `!!`
#' @importFrom tidyr pivot_longer nest unnest
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
.get_abundance_plot_data <- function(x, rank, value_select, abund_values,
                                     scale){
    if(is.null(rownames(colData(x)))){
        stop(".")
    }
    data <- assay(x,abund_values)
    if(scale){
        data <- sweep(data, 2L, colSums(data),"/")
    }
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
        unnest(cols = data) %>%
        left_join(colData(x)[value_select] %>%
                      as.data.frame %>%
                      rownames_to_column(MELT_NAME),
                  by = MELT_NAME)
    data
}

#' @importFrom tidyr pivot_wider
.order_plot_data <- function(plot_data, value_select,
                             by_abund = FALSE, decreasing = FALSE){
    if(by_abund){
        lvl <- levels(plot_data$rank)
        tmp <- plot_data %>%
            pivot_wider(id_cols = "rank", names_from = "rank",
                        values_from = "Value", values_fn = list) %>%
            as.list %>%
            lapply("[[",1L)
        o <- do.call(order,
                     c(tmp[levels(plot_data$rank)],
                       list(decreasing = decreasing)))
        lvl <- unique(plot_data[[MELT_NAME]])[o]
    } else {
        plot_data$rank <- factor(plot_data$rank, unique(plot_data$rank))
        o <- order(plot_data[[value_select]])
        lvl <- unique(plot_data[[MELT_NAME]][o])
    }
    plot_data[[MELT_NAME]] <- factor(plot_data[[MELT_NAME]], lvl)
    plot_data <- plot_data[order(plot_data[[MELT_NAME]]),]
    plot_data[[MELT_NAME]] <- as.integer(plot_data[[MELT_NAME]])
    plot_data
}

#' @importFrom ggplot2 ggplot
.get_base_plot <- function(plot_data, plot_args){
    ggplot(plot_data) +
        do.call(plot_args$theme,list())
}

#' @importFrom ggplot2 aes_string geom_bar theme element_text
#'   scale_y_continuous scale_x_continuous scale_x_discrete element_blank
.get_rel_abundance_plot_across <- function(plot, rank, across){
    aes <- aes_string(x = MELT_NAME, y = MELT_VALUES, fill = "rank")
    plot +
        geom_bar(mapping = aes,
                 stat = "identity") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0))
}

#' @importFrom ggplot2 aes_string geom_point geom_smooth theme
#'   element_blank scale_y_continuous scale_x_continuous scale_x_discrete
#'   facet_wrap
.get_rel_abundance_point_plot_across <- function(plot, rank, across){
    aes_point <- aes_string(x = MELT_NAME,
                            y = MELT_VALUES,
                            fill = "rank")
    aes_smooth <- aes_string(x = MELT_NAME,
                             y = MELT_VALUES,
                             colour = "rank")
    plot +
        geom_point(mapping = aes_point, shape = 21, colour = "grey70") +
        geom_smooth(mapping = aes_smooth, formula = y ~ x) +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        facet_wrap(. ~ rank, scales="free_y")
}

#' @importFrom ggplot2 aes_string geom_boxplot theme_bw theme
#'   element_text
.get_rel_abundance_plot_per <- function(plot, rank, per){
    aes <- aes_string(x = per, y = MELT_VALUES, colour = "rank")
    plot +
        geom_boxplot(mapping = aes) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

#' @importFrom ggplot2 guide_legend ylab xlab
.add_abundance_plot_labels <- function(plot, rank, values_select, scale,
                                       abund_values){
    if(scale){
        ylab <- paste0("Relative abundance - ",abund_values," (%)")
    } else {
        ylab <- paste0("Abundance - ",abund_values)
    }
    plot <- plot +
        guides(colour = guide_legend(title = rank)) +
        ylab(ylab) +
        xlab(values_select)
    plot
}

#' @importFrom ggplot2 theme
.add_abundance_legend <- function(plot, legend){
    if (!legend) {
        plot <- plot + theme(legend.position = "none")
    }
    plot
}

################################################################################

#' @rdname plotAbundance
#' @export
setMethod("plotAbundanceOrdered", signature = c("SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], across = NULL,
             abund_values = "relabundance", decreasing = TRUE, scale = FALSE,
             type = c("bar","point"), legend = TRUE, ...){
        # input checks
        if(!.is_non_empty_string(rank)){
          stop("'rank' must be an non empty single character value.",
               call. = FALSE)
        }
        if(nrow(x) == 0L){
          stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        if(!.is_a_bool(decreasing)){
          stop("'decreasing' must be TRUE or FALSE.", call. = FALSE)
        }
        .check_abund_values(abund_values, x)
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        across <- .norm_across(across, x)
        if(!.is_a_bool(scale)){
          stop("'scale' must be TRUE or FALSE.", call. = FALSE)
        }
        type <- match.arg(type, c("bar","point"))
        if(!.is_a_bool(legend)){
          stop("'legend' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        plot_args <- list(theme = theme_bw)
        plot_data <- .get_abundance_plot_data(x, rank, across, abund_values,
                                            scale)
        plot_data <- .order_plot_data(plot_data, across, by_abund = TRUE,
                                    decreasing = decreasing)
        #
        plot <- .get_base_plot(plot_data, plot_args)
        if(type == "point" ||
         (!scale && length(unique(plot_data[[MELT_NAME]])) == 1L)){
          plot <- .get_rel_abundance_point_plot_across(plot, rank, across)
        } else {
          plot <- .get_rel_abundance_plot_across(plot, rank, across)
        }
        plot <- .add_abundance_plot_labels(plot, rank, across, scale,
                                         abund_values)
        plot <- .add_abundance_legend(plot, legend)
        plot
    }
)
