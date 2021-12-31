#' Plot Series
#'
#' This function plots series data.
#'
#' @param object a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object.
#'
#' @param abund_values a single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   plotted. (default: \code{abund_values = "counts"})
#'
#' @param x a single character value for selecting the column from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{ColData}} that
#'   will specify values of x-axis.
#'  
#' @param y a single character value for selecting the taxa from
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{rownames}}.
#'   This parameter specifies taxa whose abundances will be plotted.
#'  
#' @param rank a single character value defining a taxonomic rank, that is used 
#'   to agglomerate the data. Must be a value of \code{taxonomicRanks()} 
#'   function.
#'  
#' @param colour_by a single character value defining a taxonomic rank, that is used to
#'   color plot. Must be a value of \code{taxonomicRanks()} function.
#' 
#' @param linetype_by a single character value defining a taxonomic rank, that
#'   is used to divide taxa to different line types. Must be a value of
#'   \code{taxonomicRanks()} function.
#' 
#' @param size_by a single character value defining a taxonomic rank, that is
#'   used to divide taxa to different line size types. Must be a value of
#'   \code{taxonomicRanks()} function.
#'   
#' @param ... additional parameters for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}
#'
#' @details
#' This function creates series plot, where x-axis includes e.g. time points, and
#' y-axis abundances of selected taxa.
#'
#' @return 
#' A \code{ggplot2} object 
#'
#' @name plotSeries
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#'
#' @examples
#' library(mia)
#' object <- microbiomeDataSets::SilvermanAGutData()
#' # Plots 2 most abudant taxa, which are colore by their family
#' plotSeries(object,
#'            x = "DAY_ORDER",
#'            y = getTopTaxa(object, 2),
#'            colour_by = "Family")
#' 
#' # Counts relative abundances
#' object <- transformCounts(object, method = "relabundance")
#' 
#' # Selects taxa
#' taxa <- c("seq_1", "seq_2", "seq_3", "seq_4", "seq_5")
#' 
#' # Plots relative abundances of phylums
#' plotSeries(object[taxa,],
#'            x = "DAY_ORDER", 
#'            colour_by = "Family",
#'            linetype_by = "Phylum",
#'            abund_values = "relabundance")
#' 
#' # In addition to 'colour_by' and 'linetype_by', 'size_by' can also be used to group taxa.
#' plotSeries(object,
#'            x = "DAY_ORDER", 
#'            y = getTopTaxa(object, 5), 
#'            colour_by = "Family",
#'            size_by = "Phylum",
#'            abund_values = "counts")
NULL

#' @rdname plotSeries
#' @export
setGeneric("plotSeries", signature = c("object"),
           function(object,
                    x,
                    y = NULL,
                    rank = NULL,
                    colour_by = NULL,
                    size_by = NULL,
                    linetype_by = NULL,
                    abund_values = "counts",
                    ...)
               standardGeneric("plotSeries"))


#' @rdname plotSeries
#' @importFrom SummarizedExperiment colData
#' @export
setMethod("plotSeries", signature = c(object = "SummarizedExperiment"),
    function(object,
             x,
             y = NULL,
             rank = NULL,
             colour_by = NULL,
             size_by = NULL,
             linetype_by = NULL,
             abund_values = "counts",
             ...){
        ###################### Input check #######################
        # Checks abund_values
        .check_assay_present(abund_values, object)
        
        # Checks X
        if( !.is_a_string(x) ||
            !(x %in% names(colData(object))) ){
            stop("'x' must be a name of column of colData(object)",
                 call. = FALSE)
        }
        
        # If rank is not null, data will be agglomerated by rank
        if( !is.null(rank) ){
            # Check rank
            .check_taxonomic_rank(rank, object)
            
            # Agglomerates the object
            object <- agglomerateByRank(object, rank = rank)
        }
        
        # Checks Y
        # If Y is not null, user has specified it
        if (!is.null(y)){
            if(!is.character(y) ||
               !all( y %in% rownames(object))){
              stop("'y' must be in rownames(x). \n If 'rank' was used, check ",
                   "that 'y' matches agglomerated data.", 
                   call. = FALSE)
            }
            # Select taxa that user has specified
            object <- object[y,]
        }
        ###################### Input check end ####################
        
        # Gets assay data
        assay <- .get_assay_data(object, abund_values)
        
        # Fetches sample and features data as a list. 
        vis_out <- .incorporate_series_vis(object,
                                           x,
                                           colour_by,
                                           linetype_by,
                                           size_by)
        series_data <- vis_out$series_data
        feature_data <- vis_out$feature_data
        x <- vis_out$x
        colour_by <- vis_out$colour_by
        linetype_by <- vis_out$linetype_by
        size_by <- vis_out$size_by
        
        # Melts the data
        plot_data <- .melt_series_data(assay,
                                       series_data,
                                       feature_data)
        xlab <- paste0(x)
        ylab <- paste0(abund_values)
        
        # Plots the data
        .series_plotter(plot_data, 
                        xlab = xlab,
                        ylab = ylab,
                        colour_by = colour_by,
                        linetype_by = linetype_by,
                        size_by = size_by,
                        ...)
    }
)

################## HELP FUNCTIONS ##########################

.get_assay_data <- function(object, abund_values){
    # Gets warning or error if too many taxa are selected. 
    if( length(rownames(object)) > 20 ){
        stop("Over 20 taxa selected. 20 or under allowed.", call. = FALSE)
    } else if ( length(rownames(object)) > 10 ){
        warning("Over 10 taxa selected.", call. = FALSE)
    }
    
    # Retrieves the assay
    assay <- assay(object, abund_values, withDimnames = TRUE)
    
    # Gets rownames
    rownames(assay) <- rownames(object)
    return(assay)
}

.incorporate_series_vis <- function(object, x, colour_by, linetype_by, size_by){
    
    # This variable is set by defaults
    series_data <- retrieveCellInfo(object, x, search = "colData")
    x <- series_data$name
    series_data <- series_data$value
    
    # This variables are optional
    row_vars <- c(colour_by = colour_by,
                  linetype_by = linetype_by,
                  size_by = size_by)
    colour_by <- NULL
    linetype_by <- NULL
    size_by <- NULL
    feature_data <- NULL
    if(!is.null(row_vars)){
        feature_data <- list()
        for(i in seq_along(row_vars)){
            # get data
            feature_data[[i]] <- retrieveFeatureInfo(object, row_vars[i],
                                                     search = "rowData")
            feature_info_name <- feature_data[[i]]$name
            # mirror back variable name, if a partial match was used
            var_name <- names(row_vars)[i]
            assign(var_name, feature_info_name)
            # rename columns by their usage
            feature_data[[i]]$name <- var_name
        }
        # squach the feature data
        if(length(feature_data) > 0L){
            names <- vapply(feature_data,"[[",character(1),"name")
            data <- lapply(feature_data,"[[","value")
            feature_data <- data.frame(data)
            colnames(feature_data) <- names
            rownames(feature_data) <- rownames(object)
        }
    }
    return(list(series_data = series_data,
                feature_data = feature_data,
                x = x,
                colour_by = colour_by,
                linetype_by = linetype_by,
                size_by = size_by))
}

#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join group_by summarize ungroup
#' @importFrom stats sd
.melt_series_data <- function(assay, series_data, feature_data){
    colnames(assay) <- seq_len(ncol(assay))
    # Melts assay table
    melted_data <- as.data.frame(assay) %>% 
        rownames_to_column("feature") %>% 
        pivot_longer(-c("feature"),
                     names_to = "sample",
                     values_to = "Y")
    # joing the series data
    melted_data <- melted_data %>%
        left_join(data.frame(sample = colnames(assay),
                             X = series_data),
                  by = "sample") %>%
        select(!sym("sample"))
    # if replicates are present calculate sd
    if(anyDuplicated(melted_data$X)){
        melted_data <- melted_data %>% 
            group_by(!!sym("X"),!!sym("feature")) %>%
            summarize(sd = sd(.data$Y, na.rm = TRUE),
                      Y = mean(.data$Y, na.rm = TRUE)) %>%
            ungroup()
    }
    # join the feature data
    if(!is.null(feature_data)){
        feature_data <- feature_data %>%
            rownames_to_column("feature")
        melted_data <- melted_data %>%
            left_join(feature_data,
                      by = "feature")
    }
    melted_data
}

.series_plotter <- function(plot_data,
                            xlab = NULL,
                            ylab = NULL,
                            colour_by = NULL,
                            linetype_by = NULL,
                            size_by = NULL,
                            add_legend = TRUE,
                            line_alpha = 1,
                            line_type = NULL,
                            line_width = 1,
                            line_width_range = c(0.5,3),
                            ribbon_alpha = 0.3){
    # fall back for feature grouping
    if(is.null(colour_by)){
        colour_by <- "feature"
        plot_data$colour_by <-  plot_data$feature
    }
    # Creates a "draft" of a plot
    plot_out <- ggplot(plot_data,
                       aes_string(x = "X", y = "Y")) +
        labs(x = xlab, y = ylab)
    # if sd column is present add a ribbon
    if(!is.null(plot_data$sd)){
        ribbon_args <- .get_ribbon_args(colour_by = colour_by,
                                        alpha = ribbon_alpha)
        plot_out <- plot_out +
            do.call(geom_ribbon, ribbon_args$args)
    }
    # Fetches arguments fpr geom_line
    line_args <- .get_line_args(colour_by = colour_by,
                                linetype_by = linetype_by,
                                size_by = size_by,
                                alpha = line_alpha,
                                linetype = line_type,
                                size = line_width)
    # Adds arguments to the plot
    plot_out <- plot_out +
        do.call(geom_line, line_args$args)
    # apply line_width_range
    if (!is.null(size_by)) {
        if(is.numeric(plot_data$size_by)){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out +
            SIZEFUN(range = line_width_range)
    }
    # Resolves the colours
    plot_out <- .resolve_plot_colours(plot_out,
                                      plot_data$colour_by,
                                      colour_by,
                                      fill = FALSE)
    if(!is.null(plot_data$sd)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          plot_data$colour_by,
                                          colour_by,
                                          fill = TRUE)
    } 
    # add additional guides
    plot_out <- .add_extra_line_guide(plot_out, linetype_by, size_by)
    # To choose if legend is kept, and its position
    plot_out <- .add_legend(plot_out, add_legend)
    plot_out <- plot_out +
        theme_classic()
    plot_out
}

.add_extra_line_guide <- function(plot_out, linetype_by, size_by) {
    guide_args <- list()
    if (!is.null(linetype_by)) {
        guide_args$linetype <- guide_legend(title = linetype_by)
    }
    if (!is.null(size_by)) {
        guide_args$size <- guide_legend(title = size_by)
    }
    if (length(guide_args)) {
        plot_out <- plot_out + do.call(guides, guide_args)
    }
    return(plot_out)
}
