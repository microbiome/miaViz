#' Plotting abundance data
#'
#' \code{plotAbundance} plots the abundance on a selected taxonomic rank.
#' Since this probably makes sense only for relative abundance data, the
#' assay used by default is expected to be in the slot \sQuote{relabundance}.
#' If only \sQuote{counts} is present, the relative abundance is computed.
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
#' @param features a single \code{character} value defining a column from 
#'   \code{colData} to be plotted below the abundance plot.
#'   Continuous numeric values will be plotted as point, whereas factors and
#'   character will be plotted as colour-code bar. (default: \code{features =
#'   NULL})
#'   
#' @param order_rank_by How to order abundance value: By name (\dQuote{name}) 
#' for sorting the taxonomic labels alphabetically, by abundance (\dQuote{abund}) to
#' sort by abundance values or by a reverse order of abundance values (\dQuote{revabund}).
#'  
#'   
#' @param order_sample_by A single character value from the chosen rank of abundance
#'   data or from \code{colData} to select values to order the abundance
#'   plot by. (default: \code{order_sample_by = NULL})
#'   
#' @param decreasing TRUE or FALSE: If the \code{order_sample_by} is defined and the
#'   values are numeric, should the values used to order in decreasing or
#'   increasing fashion? (default: \code{decreasing = FALSE})
#'   
#' @param use_relative \code{TRUE} or \code{FALSE}: Should the relative values
#'   be calculated? (default: \code{use_relative = TRUE})
#'
#' @param layout Either \dQuote{bar} or \dQuote{point}. 
#' 
#' @param one_facet Should the plot be returned in on facet or split into 
#'   different facet, one facet per different value detect in \code{rank}. If
#'   \code{features} or \code{order_sample_by} is not \code{NULL}, this setting will
#'   be disregarded.
#' 
#' @param ncol,scales if \code{one_facet = FALSE}, \code{ncol} defines many 
#'   columns should be for plotting the different facets and \code{scales} is
#'   used to define the behavior of the scales of each facet. Both values are 
#'   passed onto \code{\link[ggplot2:facet_wrap]{facet_wrap}}.
#' 
#' @param ... additional parameters for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}
#'
#' @return 
#' a \code{\link[ggplot2:ggplot]{ggplot}} object or list of two 
#' \code{\link[ggplot2:ggplot]{ggplot}} objects, if `features` are added to 
#' the plot. 
#'
#' @name plotAbundance
#'
#' @examples
#' data(GlobalPatterns, package="mia")
#' se <- GlobalPatterns
#' 
#' ## Plotting abundance using the first taxonomic rank as default
#' plotAbundance(se, abund_values="counts")
#' 
#' ## Using "Phylum" as rank
#' plotAbundance(se, abund_values="counts", rank = "Phylum", add_legend = FALSE)
#' 
#' ## If rank is set to NULL plotAbundance behaves like plotExpression
#' plotAbundance(se, abund_values="counts", rank = NULL,
#'            features = head(rownames(se)))
#'   
#' ## A feature from colData or taxon from chosen rank can be used for ordering samples.
#' plotAbundance(se, abund_values="counts", rank = "Phylum",
#'            order_sample_by = "Bacteroidetes")
#' 
#' ## Features from colData can be plotted together with abundance plot.
#' # Returned object is a list that includes two plot; other visualizes abundance
#' # other features. 
#' plot <- plotAbundance(se, abund_values = "counts", rank = "Phylum",
#'                    features = "SampleType")
#' \donttest{
#' # These two plots can be combined with wrap_plots function from patchwork package
#' library(patchwork)
#' wrap_plots(plot, ncol = 1)
#' }
#' 
#' ## Compositional barplot with top 5 taxa and samples sorted by "Bacteroidetes"
#' 
#' # Getting top taxa on a Phylum level
#' se <- relAbundanceCounts(se)
#' se_phylum <- agglomerateByRank(se, rank ="Phylum", onRankOnly=TRUE)
#' top_taxa <- getTopTaxa(se_phylum,top = 5, abund_values = "relabundance")
#' 
#' # Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
#' phylum_renamed <- lapply(rowData(se)$Phylum,
#'                        function(x){if (x %in% top_taxa) {x} else {"Other"}})
#' rowData(se)$Phylum <- as.character(phylum_renamed)
#' 
#' # Compositional barplot
#' plotAbundance(se, abund_values="relabundance", rank = "Phylum",
#'            order_rank_by="abund", order_sample_by = "Bacteroidetes")
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
            features = NULL,
            order_rank_by = c("name","abund","revabund"),
            order_sample_by = NULL,
            decreasing = TRUE,
            use_relative = TRUE,
            layout = c("bar","point"),
            one_facet = TRUE,
            ncol = 2,
            scales = "fixed",
            abund_values = "counts",
            ...){
        # input checks
        if(nrow(x) == 0L){
            stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        .check_assay_present(abund_values, x)
        # if rank is set to NULL, default to plotExpression 
        if(is.null(rank)){
            plot <- plotExpression(x, features = features, 
                                exprs_values = abund_values,
                                one_facet = one_facet,
                                ncol = ncol, scales = scales, ...)
            ylab <- gsub("Expression","Abundance",plot$labels$y)
            plot <- plot +
                ylab(ylab)
            return(plot)
        }
        ############################# INPUT CHECK #############################
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be an non empty single character value.",
                call. = FALSE)
        }
        if(!.is_a_bool(use_relative)){
            stop("'use_relative' must be TRUE or FALSE.",
                call. = FALSE)
        }
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        layout <- match.arg(layout, c("bar","point"))
        order_rank_by <- match.arg(order_rank_by, c("name","abund","revabund"))
        .check_abund_plot_args(one_facet = one_facet,
                            ncol = ncol)
        if( !is.null(features) ){
            features <- match.arg(features, colnames(colData(x)))
        }
        ########################### INPUT CHECK END ###########################
        abund_data <- .get_abundance_data(x, rank, abund_values, order_rank_by,
                                        use_relative)
        order_sample_by <- .norm_order_sample_by(order_sample_by,
                                                unique(abund_data$colour_by),
                                                x)
        features_data <- NULL
        if(!is.null(features) || !is.null(order_sample_by)){
            features_data <- .get_features_data(features, order_sample_by, x)
        }
        if(!is.null(order_sample_by)){
            order_out <- .order_abund_feature_data(abund_data, features_data,
                                                order_sample_by, decreasing)
            abund_data <- order_out$abund_data
            features_data <- order_out$features_data
            rm(order_out)
        }
        plot_out <- .abund_plotter(abund_data,
                                xlab = "Samples",
                                ylab = "Rel. Abundance",
                                colour_by = rank,
                                layout = layout,
                                ...)
        if(!is.null(features_data)){
            plot_feature_out <- .features_plotter(features_data,
                                                order_sample_by,
                                                xlab = "Samples",
                                                ...)
            plot_out <- c(list(abundance = plot_out), plot_feature_out)
        } else {
            if (!one_facet) {
                plot_out <- plot_out + 
                    facet_wrap(~colour_by, ncol = ncol, scales = scales)
            }
        }
        # Checks if the list is a ggplot object or regular list of ggplot objects
        if( !is.ggplot(plot_out) ){
            # If features is specified, then only abundance and features plots are 
            # returned as a list. If it is not, then only abundance plot is returned.
            if( !is.null(features) ){
                plot_out <- list(abundance = plot_out[["abundance"]], plot_out[[features]])
                # Assigns the names back
                names(plot_out) <- c("abundance", features)
            } else{
                plot_out <- plot_out[["abundance"]]
            }
        }
        plot_out
    }
)


MELT_NAME <- "Sample"
MELT_VALUES <- "Value"

#' @importFrom SummarizedExperiment rowData assay
#' @importFrom dplyr mutate group_by summarize rename
#' @importFrom tidyr pivot_longer nest unnest
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
.get_abundance_data <- function(x, rank, abund_values, order_rank_by = "name",
                                use_relative = TRUE){
    data <- assay(x, abund_values, withDimnames = TRUE)
    if(use_relative){
        data <- mia:::.calc_rel_abund(data)
    }
    if(is.null(colnames(x))){
        colnames(data) <- paste0("Sample",seq_len(ncol(x)))
    }
    merge_FUN <- function(data){
        data %>%
            group_by(!!sym(MELT_NAME)) %>%
            summarize(!!sym(MELT_VALUES) := sum(!!sym(MELT_VALUES)))
    }
    # enable conversion to data.frame for non-matrix assays, e.g. sparseMatrices
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    data <- data %>%
        as.data.frame() %>%
        mutate(rank = factor(rowData(x)[,rank], unique(rowData(x)[,rank]))) %>%
        pivot_longer(cols = !.data$rank,
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
    # order values
    if(order_rank_by == "name"){
        lvl <- levels(data$colour_by)
        lvl <- lvl[order(lvl)]
    } else if(order_rank_by %in% c("abund","revabund")){
        o <- data %>% 
            select(!.data$X) %>% 
            group_by(.data$colour_by) %>% 
            summarize(sum = sum(.data$Y))
        decreasing <- ifelse(order_rank_by == "abund",TRUE,FALSE)
        lvl <- o[order(o$sum, decreasing = decreasing),] %>% 
            pull(.data$colour_by) %>%
            as.character()
    } else {
        stop(".")
    }
    data$colour_by <- factor(data$colour_by, lvl)
    data <- data[order(data$colour_by),]
    #
    data
}

.norm_order_sample_by <- function(order_sample_by, factors, x){
    if(is.null(order_sample_by)){
        return(order_sample_by)
    }
    msg <- paste0("'order_sample_by' must be a single non-empty character value, ",
                "either present in the abundance data as group variable or ",
                "in the column data of 'x'. (The abundance data takes ",
                "precedence)")
    if(!.is_non_empty_string(order_sample_by)){
        stop(msg, call. = FALSE)
    }
    if(!(order_sample_by %in% factors)){
        tmp <- try({retrieveCellInfo(x, order_sample_by, search = "colData")},
                silent = TRUE)
        if(is(tmp,"try-error")){
            stop(msg, call. = FALSE)
        } else {
            order_sample_by <- tmp$name
        }
    }
    order_sample_by
}

.get_feature_data <- function(x, by){
    tmp <- try({retrieveCellInfo(x, by, search = "colData")}, silent = TRUE)
    if(is(tmp,"try-error")){
        tmp <- NULL
    }
    tmp
}

.get_features_data <- function(features, order_sample_by, x){
    features <- unique(c(order_sample_by,features))
    features_data <- lapply(features,
                            .get_feature_data,
                            x = x)
    non_empty <- !vapply(features_data, is.null, logical(1))
    features_data <- features_data[non_empty]
    if(length(features_data) == 0){
        return(NULL)
    }
    values <- lapply(features_data, "[[", "value")
    names <- lapply(features_data, "[[", "name")
    features_data <- data.frame(values)
    colnames(features_data) <- names
    if(!is.null(colnames(x))){
        rownames(features_data) <- colnames(x)
    } else {
        rownames(features_data) <- paste0("Sample",seq_len(ncol(x)))
    }
    features_data
}

#' @importFrom dplyr pull
.order_abund_feature_data <- function(abund_data, features_data,
                                    order_sample_by, decreasing = TRUE){
    if(!is.null(order_sample_by)){
        lvl <- levels(abund_data$X)
        if(is.null(features_data) || 
            !(order_sample_by %in% colnames(features_data))){
            # presort by rank value
            lvl_tmp <- levels(abund_data$colour_by)
            lvl_tmp <- c(order_sample_by, lvl_tmp[!(lvl_tmp %in% order_sample_by)])
            abund_data$colour_by <- factor(abund_data$colour_by, lvl_tmp)
            #
            data <- abund_data[abund_data$colour_by %in% order_sample_by,] %>%
                pull("Y")
        } else {
            data <- features_data[,order_sample_by]
        }
        if(is.factor(data)){
            o <- order(data, decreasing = !decreasing)
        } else {
            o <- order(data, decreasing = decreasing)
        }
        lvl <- lvl[o]
        # reset lvls and reorder
        abund_data$X <- factor(abund_data$X, lvl)
        abund_data <- abund_data[order(abund_data$colour_by, abund_data$X),]
        if(!is.null(features_data)){
            o <- order(factor(rownames(features_data), lvl))
            # If features and order_sample_by are the same, there is only one column.
            # One column is converted to vector which is why it is converted back
            # to data frame which is expected in next steps.
            if(ncol(features_data) == 1) {
                colname <- colnames(features_data)
                features_data <- as.data.frame(features_data[o,])
                colnames(features_data) <- colname
            } else{
                features_data <- features_data[o,]
            }
        }
    }
    list(abund_data = abund_data, features_data = features_data)
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
                        add_border = NULL,
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
                                alpha = bar_alpha,
                                add_border = add_border,
                                n = length(unique(object$X)))
        plot_out <- plot_out +
            do.call(geom_bar, c(abund_out$args, list(stat="identity"))) +
            scale_y_continuous(expand = c(0,0))

    } else {
        abund_out <- .get_point_args(colour_by,
                                    shape_by = NULL,
                                    size_by = NULL,
                                    alpha = point_alpha,
                                    size = point_size)
        abund_out$border <- TRUE
        plot_out <- plot_out +
            do.call(geom_point, abund_out$args)
    }
    # adjust point colours
    if(!is.null(colour_by)){
        if(abund_out$border){
            # resolve the colour for the line colours
            plot_out <- .resolve_plot_colours(plot_out,
                                            object$colour_by,
                                            colour_by,
                                            fill = FALSE)
        }
        # resolve the color for fill
        plot_out <- .resolve_plot_colours(plot_out,
                                        object$colour_by,
                                        colour_by,
                                        fill = TRUE)
    }
    plot_out <- plot_out +
        theme_classic()
    # add legend
    plot_out <- .add_legend(plot_out, add_legend)
    # flip
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    plot_out
}

#' @importFrom ggplot2 ggplot aes_string labs geom_point geom_raster
.feature_plotter <- function(feature_data,
                            name,
                            xlab,
                            flipped,
                            add_legend,
                            add_x_text,
                            point_alpha,
                            point_size){
    if(is.factor(feature_data$Y)){
        feature_data$colour_by <- feature_data$Y
        feature_data$Y <- ""
        colour_by <- unique(feature_data$feature_name)
    }
    feature_plot_out <- ggplot(feature_data, aes_string(x="X", y="Y")) +
        labs(x = xlab, y = name)
    if(length(unique(feature_data$Y)) == 1L){
        feature_out <- .get_bar_args(colour_by,
                                    alpha = point_alpha,
                                    add_border = FALSE)
        feature_plot_out <- feature_plot_out +
            do.call(geom_raster, feature_out$args) +
            scale_y_discrete(expand = c(0,0))
        feature_plot_out <- .resolve_plot_colours(feature_plot_out,
                                                feature_data$colour_by,
                                                colour_by,
                                                fill = TRUE)
        legend_pos <- "bottom"
    } else {
        feature_out <- .get_point_args(NULL,
                                    shape_by = NULL,
                                    size_by = NULL,
                                    alpha = point_alpha,
                                    size = point_size)
        feature_plot_out <- feature_plot_out +
            do.call(geom_point, feature_out$args)
        legend_pos <- "right"
    }
    feature_plot_out <- feature_plot_out +
        theme_classic()
    # add legend
    feature_plot_out <- .add_legend(feature_plot_out, add_legend, legend_pos)
    # flip
    feature_plot_out <- .flip_plot(feature_plot_out, flipped, add_x_text)
    feature_plot_out
}

.features_plotter <- function(features_data,
                            order_sample_by = NULL,
                            xlab = NULL,
                            flipped = FALSE,
                            add_legend = TRUE,
                            add_x_text = FALSE,
                            point_alpha = 1,
                            point_size = 2,
                            ...){
    names <- colnames(features_data)
    features_data <- lapply(names, 
                            function(col){
                                data.frame(X = factor(rownames(features_data),
                                                    rownames(features_data)),
                                                    feature_name = col,
                                                    Y = features_data[,col])
                            })
    names(features_data) <- names
    plots_out <- mapply(.feature_plotter,
                        features_data,
                        names(features_data),
                        MoreArgs = list(xlab = xlab,
                                        flipped = flipped,
                                        add_legend = add_legend,
                                        add_x_text = add_x_text,
                                        point_alpha = point_alpha,
                                        point_size = point_size),
                        SIMPLIFY = FALSE)
    names(plots_out) <-  names(features_data)
    if(!is.null(order_sample_by)){
        reorder <- c(order_sample_by,
                    names(plots_out)[!(names(plots_out) %in% order_sample_by)])
        m <- match(reorder,names(plots_out))
        m <- m[!is.na(m)]
        plots_out <- plots_out[m]
    }
    plots_out
}
