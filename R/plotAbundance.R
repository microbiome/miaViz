#' Plotting abundance data
#'
#' \code{plotAbundance} plots the abundance on a selected column in rowData.
#' Since this probably makes sense only for relative abundance data, the
#' assay used by default is expected to be in the slot \sQuote{relabundance}.
#' If only \sQuote{counts} is present, the relative abundance is computed.
#'
#' Subsetting to rows of interested and ordering of those is expected to be done
#' outside of this functions, e.g. \code{x[1:2,]}. This will plot data of all
#' col.var present.
#'
#' @param x a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object.
#'
#' @param group \code{Character scalar}. Defines a group to use. Must be a value 
#' of \code{colnames(rowData(x))}. (Default: \code{NULL})
#'   
#' @param rank Deprecated. Use \code{group} instead.
#'
#' @param assay.type \code{Character scalar} value defining which assay data to
#'   use. (Default: \code{"relabundance"})
#'   
#' @param assay_name Deprecate. Use \code{assay.type} instead.
#'   
#' @param col.var \code{Character scalar}. Selects a column from 
#'   \code{colData} to be plotted below the abundance plot.
#'   Continuous numeric values will be plotted as point, whereas factors and
#'   character will be plotted as colour-code bar. (Default: \code{NULL})
#'   
#' @param features Deprecated. Use \code{col.var} instead.
#'   
#' @param order.row.by \code{Character scalar}. How to order abundance value: By name (\dQuote{name}) 
#' for sorting the taxonomic labels alphabetically, by abundance (\dQuote{abund}) to
#' sort by abundance values or by a reverse order of abundance values (\dQuote{revabund}).
#' 
#' @param order_rank_by Deprecated. Use \code{order.row.by} instead.  
#'   
#' @param order.col.by \code{Character scalar}. from the chosen rank of abundance
#'   data or from \code{colData} to select values to order the abundance
#'   plot by. (Default: \code{NULL})
#'   
#' @param order_sample_by Deprecated. Use \code{order.col.by} instead.
#'   
#' @param decreasing \code{Logical scalar}. If the \code{order.col.by} is defined and the
#'   values are numeric, should the values used to order in decreasing or
#'   increasing fashion? (Default: \code{FALSE})
#'
#' @param layout \code{Character scalar}. Either \dQuote{bar} or \dQuote{point}. 
#' 
#' @param one.facet \code{Logical scalar}. Should the plot be returned in on facet or split into 
#'   different facet, one facet per different value detect in \code{rank}. If
#'   \code{col.var} or \code{order.col.by} is not \code{NULL}, this setting will
#'   be disregarded. (Default: \code{TRUE})
#'   
#' @param one_facet Deprecated. Use \code{one.facet} instead.
#' 
#' @param ncol \code{Numeric scalar}. if \code{one.facet = FALSE}, \code{ncol} defines many 
#'   columns should be for plotting the different facets. (Default: \code{2})
#'   
#' @param scales \code{Character scalar}. Defines the behavior of the scales of each facet. Both values are 
#'   passed onto \code{\link[ggplot2:facet_wrap]{facet_wrap}}. (Default: \code{"fixed"})
#' 
#' @param ... additional parameters for plotting.
#'   \itemize{
#'   \item \code{as.relative} \code{Character scalar}. Should the relative values
#'   be calculated? (Default: \code{FALSE})
#' }
#' See \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}
#'
#' @return 
#' a \code{\link[ggplot2:ggplot]{ggplot}} object or list of two 
#' \code{\link[ggplot2:ggplot]{ggplot}} objects, if `col.var` are added to 
#' the plot. 
#'
#' @name plotAbundance
#'
#' @examples
#' data(GlobalPatterns, package="mia")
#' tse <- GlobalPatterns
#' 
#' ## If rank is set to NULL (default), agglomeration is not done. However, note
#' ## that there is maximum number of rows that can be plotted. That is why
#' ## we take sample from the data.
#' set.seed(26348)
#' sample <- sample(rownames(tse), 20)
#' tse_sub <- tse[sample, ]
#' # Apply relative transformation
#' tse_sub <- transformAssay(tse_sub, method = "relabundance")
#' plotAbundance(tse_sub, assay.type = "relabundance")
#' 
#' ## Plotting counts using the first taxonomic rank as default
#' plotAbundance(
#'     tse, assay.type="counts", rank = "Phylum") +
#'     labs(y="Counts")
#' 
#' ## Using "Phylum" as rank. Apply relative transformation to "counts" assay.
#' plotAbundance(
#'     tse, assay.type="counts", rank = "Phylum", add_legend = FALSE,
#'     as.relative = TRUE)
#' 
#' # Apply relative transform
#' tse <- transformAssay(tse, method = "relabundance")
#'   
#' ## A feature from colData or taxon from chosen rank can be used for ordering
#' ## samples.
#' plotAbundance(tse, assay.type="relabundance", rank = "Phylum",
#'            order.col.by = "Bacteroidetes")
#' 
#' ## col.var from colData can be plotted together with abundance plot.
#' # Returned object is a list that includes two plot; other visualizes
#' ## abundance other col.var. 
#' plot <- plotAbundance(tse, assay.type = "relabundance", rank = "Phylum",
#'                    col.var = "SampleType")
#' \donttest{
#' # These two plots can be combined with wrap_plots function from patchwork
#' # package
#' library(patchwork)
#' wrap_plots(plot, ncol = 1)
#' }
#' 
#' ## Same plot as above but showing sample IDs as labels for the x axis on the
#' ## top plot
#' plot[[1]] <- plotAbundance(tse, assay.type = "relabundance", rank = "Phylum",
#'                            col.var = "SampleType", add.legend = FALSE,
#'                            add.x.text = TRUE)[[1]] +
#'                            theme(axis.text.x = element_text(angle = 90)) 
#' \donttest{
#' wrap_plots(plot, ncol = 1, heights = c(0.8,0.2))
#' }
#' 
#' ## Compositional barplot with top 5 taxa and samples sorted by
#' ## "Bacteroidetes"
#' 
#' # Getting top taxa on a Phylum level
#' tse <- transformAssay(tse, method="relabundance")
#' tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
#' top_taxa <- getTop(tse_phylum,top = 5, assay.type = "relabundance")
#' 
#' # Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
#' phylum_renamed <- lapply(rowData(tse)$Phylum,
#'                        function(x){if (x %in% top_taxa) {x} else {"Other"}})
#' rowData(tse)$Phylum <- as.character(phylum_renamed)
#' 
#' # Compositional barplot
#' plotAbundance(tse, assay.type="relabundance", rank = "Phylum",
#'            order.row.by="abund", order.col.by = "Bacteroidetes")
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
            group = rank,
            rank = NULL,
            col.var = features,
            features = NULL,
            order.row.by = order_rank_by,
            order_rank_by = c("name","abund","revabund"),
            order.col.by = order_sample_by,
            order_sample_by = NULL,
            decreasing = TRUE,
            layout = c("bar","point"),
            one.facet = one_facet,
            one_facet = TRUE,
            ncol = 2,
            scales = "fixed",
            assay.type = assay_name, assay_name = "counts",
            ...){
        ############################# INPUT CHECK #############################
        if(nrow(x) == 0L){
            stop("No data to plot. nrow(x) == 0L.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        if(!.is_non_empty_string(group) && !is.null(group)){
            stop("'group' must be an non empty single character value or NULL.",
                call. = FALSE)
        }
        if(!is.null(group)){
            if(!(group %in% colnames(rowData(x)))){
                stop("'group' must be a column from rowData .",
                     call. = FALSE)
            }
        }
        .check_for_taxonomic_data_order(x)
        layout <- match.arg(layout, c("bar","point"))
        order.row.by <- match.arg(order.row.by, c("name","abund","revabund"))
        .check_abund_plot_args(one_facet = one.facet,
                            ncol = ncol)
        if( !is.null(col.var) ){
            col.var <- match.arg(col.var, colnames(colData(x)))
        }
        ########################### INPUT CHECK END ###########################
        # Get the abundance data to be plotted. Agglomerate and apply relative
        # transformation if specified.
        abund_data <- .get_abundance_data(
            x, group, assay.type, order.row.by, ...)
        # If group was NULL, then the data was not agglomerated. The group is
        # still used in coloring (passed to colour_by parameter in
        # .abund_plotter), which is why we adjust the value of it to apply
        # coloring in (NULL means that coloring is not applied).
        group <- ifelse(is.null(group), "Feature", group)
        # Order columns
        order_col_by <- .norm_order_sample_by(
            order.col.by, unique(abund_data$colour_by), x)
        # Get additional column metadata to be plotted
        features_data <- NULL
        if(!is.null(col.var) || !is.null(order_col_by)){
            features_data <- .get_features_data(col.var, order_col_by, x)
        }
        # Order the whole data to follow user specified ordering
        if(!is.null(order_col_by)){
            order_out <- .order_abund_feature_data(
                abund_data, features_data, order_col_by, decreasing)
            abund_data <- order_out$abund_data
            features_data <- order_out$features_data
        }
        # Create the main plot
        plot_out <- .abund_plotter(abund_data,
                                colour_by = group,
                                layout = layout,
                                ...)
        # Create the column metadata plot and create a list from plots
        if(!is.null(features_data)){
            plot_feature_out <- .features_plotter(
                features_data, order.col.by, ...)
            plot_out <- c(list(abundance = plot_out), plot_feature_out)
        } else {
            # Whether to split the main plot to multiple facets. This is
            # disabled if user wants to plot also column metadata.
            if (!one.facet) {
                plot_out <- plot_out + 
                    facet_wrap(~colour_by, ncol = ncol, scales = scales)
            }
        }
        # Checks if the list is a ggplot object or regular list of ggplot
        # objects
        if( !is.ggplot(plot_out) ){
            # If features is specified, then only abundance and features plots
            # are returned as a list. If it is not, then only abundance plot is
            # returned.
            if( !is.null(col.var) ){
                plot_out <- list(
                    abundance = plot_out[["abundance"]], plot_out[[col.var]])
                # Assigns the names back
                names(plot_out) <- c("abundance", col.var)
            } else{
                plot_out <- plot_out[["abundance"]]
            }
        }
        return(plot_out)
    }
)

#' @importFrom dplyr group_by summarize rename
#' @importFrom mia meltSE
.get_abundance_data <- function(
        x, group, assay.type, order_rank_by = "name", as.relative = use_relative,
        use_relative = FALSE, ...){
    # Input check
    if(!.is_a_bool(as.relative)){
        stop("'as.relative' must be TRUE or FALSE.",
             call. = FALSE)
    }
    #
    # Agglomerate data if user has specified
    if (!is.null(group) && group %in% taxonomyRanks(x)) {
        x <- agglomerateByRank(x, group, ...)
        # or factor that is specified by user
    } else if (!is.null(group)) {
        x <- agglomerateByVariable(x, by = "rows", f = group, ...)
    }
    # At this point, we can check how many rows there are to plot. In practice,
    # there is a limit how many rows we can plot. If there are too many, it is
    # impossible to read the plot. Moreover, the plotting takes excessive
    # amount of time. The good limit might be somewhere around 50, but it
    # might be better to have higher maximum limit so we do not limit too much.
    max_num <- 500
    if( nrow(x) > max_num ){
        stop("The data contains more than ", max_num, " rows. The abundance ",
            "plot cannot be created. Consider subsetting/agglomeration. ",
            "(Check 'group' parameter)", call. = FALSE)
    }
    # If user wants to calculate relative abundances, apply relative transform
    # and use relative assay instead of the original assay in plotting.
    if( as.relative ){
        temp_name <- "temporary_relative_abundance"
        x <- transformAssay(
            x, assay.type = assay.type, method = "relabundance",
            name = temp_name)
        assay.type <- temp_name
    }
    # Samples must have names. In theory, TreeSE can include columns without
    # names. If that is the case, add names.
    if( is.null(colnames(x)) ){
        colnames(x) <- paste0("Sample", seq_len(ncol(x)))
    }
    # Melt TreeSE
    data <- meltSE(
        x, assay.type = assay.type, row.name = "colour_by", col.name = "X")
    # Add correct column name for abundance values
    colnames(data)[ colnames(data) == assay.type ] <- "Y"
    # Reorder so that the order follows the order of sample names. The order is
    # currently alphabetical.
    data$X <- factor(data$X, levels = colnames(x))
    # Order values
    if( order_rank_by == "name" ){
        # By default, factors follow alphabetical order. Order values, based
        # on names i.e. alphabetical order.
        lvl <- levels(data$colour_by)
        lvl <- lvl[order(lvl)]
    } else if( order_rank_by %in% c("abund","revabund") ){
        # Control the order
        decreasing <- ifelse(order_rank_by == "abund",TRUE,FALSE)
        # Get values for each feature and sum them up
        o <- data %>% 
            select(!.data$X) %>%
            group_by(.data$colour_by) %>% 
            summarize(sum = sum(.data$Y))
        # Order the data based on abundance. By default, feature that has
        # highest library size comes first.
        lvl <- o[order(o$sum, decreasing = decreasing), ] %>% 
            pull(.data$colour_by) %>%
            as.character()
    } else {
        stop(".")
    }
    # Apply the order
    data$colour_by <- factor(data$colour_by, lvl)
    data <- data[order(data$colour_by),]
    
    return(data)
}

.norm_order_sample_by <- function(order_sample_by, factors, x){
    # If user did not specify the ordering, do nnothing
    if(is.null(order_sample_by)){
        return(order_sample_by)
    }
    # Chec that the parameter is string
    msg <- paste0("'order.col.by' must be a single non-empty character value, ",
                "either present in the abundance data as group variable or ",
                "in the column data of 'x'. (The abundance data takes ",
                "precedence)")
    if(!.is_non_empty_string(order_sample_by)){
        stop(msg, call. = FALSE)
    }
    #
    # If the order variable is not found from the coloring values (samples can
    # also be order based on abundance of certain taxon), check that the
    # variable can be found from colData.
    if(!(order_sample_by %in% factors)){
        tmp <- .get_feature_data(x, order_sample_by, msg = msg)
        order_sample_by <- tmp$name
    }
    return(order_sample_by)
}

# This funtion retrives data from colData without failing (if the variable
# is not found)
.get_feature_data <- function(x, by, msg = NULL){
    # Get variable without failing
    tmp <- try({retrieveCellInfo(x, by, search = "colData")}, silent = TRUE)
    # Check if fetching was failed
    if(is(tmp,"try-error")){
        # If msg is not NULL, fail if variable is not found.
        # Otherwise, give NULL.
        if( !is.null(msg) ){
            stop(msg, call. = FALSE)
        } else{
            tmp <- NULL
        }
    }
   return(tmp)
}

.get_features_data <- function(features, order_sample_by, x){
    # Get variable (this functions fetches data for both odering and feature
    # plotting)
    features <- unique(c(order_sample_by,features))
    # Get values
    features_data <- lapply(
        features, .get_feature_data, x = x)
    # Get those values that are found
    non_empty <- !vapply(features_data, is.null, logical(1))
    features_data <- features_data[non_empty]
    if(length(features_data) == 0){
        return(NULL)
    }
    # Get values and names and create data.frame from them
    values <- lapply(features_data, "[[", "value")
    names <- lapply(features_data, "[[", "name")
    features_data <- data.frame(values)
    colnames(features_data) <- names
    # Add sample names to rows
    if(!is.null(colnames(x))){
        rownames(features_data) <- colnames(x)
    } else {
        rownames(features_data) <- paste0("Sample",seq_len(ncol(x)))
    }
    return(features_data)
}

#' @importFrom dplyr pull
.order_abund_feature_data <- function(
        abund_data, features_data, order_sample_by, decreasing = TRUE){
    # If ordering was specified
    if(!is.null(order_sample_by)){
        # Get levels
        lvl <- levels(abund_data$X)
        # If user specified taxan to be used to order the data
        if(is.null(features_data) || 
                !(order_sample_by %in% colnames(features_data))){
            # Presort by taxon value
            lvl_tmp <- levels(abund_data$colour_by)
            lvl_tmp <- c(
                order_sample_by, lvl_tmp[!(lvl_tmp %in% order_sample_by)])
            abund_data$colour_by <- factor(abund_data$colour_by, lvl_tmp)
            # Get the abundance values from certain taxon
            data <- abund_data[abund_data$colour_by %in% order_sample_by, ]
            data <- data[["Y"]]
        } else {
            # Otherwise get the order from the variable
            data <- features_data[[order_sample_by]]
        }
        # If the ordering value is factor, order based on alphabetical order.
        # Order numeric values in incerasing order.
        if(is.factor(data)){
            o <- order(data, decreasing = !decreasing)
        } else {
            o <- order(data, decreasing = decreasing)
        }
        # Reset lvls and reorder the data based on 
        lvl <- lvl[o]
        abund_data$X <- factor(abund_data$X, lvl)
        abund_data <- abund_data[order(abund_data$colour_by, abund_data$X),]
        # If the data includes also sample metadata to be plotted, reorder
        # it also.
        if(!is.null(features_data)){
            o <- order(factor(rownames(features_data), lvl))
            features_data <- features_data[o, , drop = FALSE]
        }
    }
    res <- list(abund_data = abund_data, features_data = features_data)
    return(res)
}


################################################################################
# abundance plotter


#' @importFrom ggplot2 ggplot theme_classic geom_point geom_bar coord_flip
#'   scale_y_continuous
.abund_plotter <- function(object,
        xlab = "Samples",
        ylab = paste0(ifelse(as.relative, "Rel. ", ""),"Abundance"),
        colour_by = NULL,
        layout = "bar",
        flipped = FALSE,
        add_legend = TRUE,
        add_x_text = add.x.text,
        add.x.text = FALSE,
        add_border = add.border,
        add.border = NULL,
        bar_alpha = bar.alpha,
        bar.alpha = 0.65,
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        as.relative = use_relative,
        use_relative = FALSE,
        ...
        ){
    # start plotting
    plot_out <- ggplot(object, aes(x=.data[["X"]], y=.data[["Y"]])) +
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
    # Adjust theme
    plot_out <- plot_out +
        theme_classic()
    # Remove legend if speicified
    plot_out <- .add_legend(plot_out, add_legend)
    # Flip the plot if specified
    plot_out <- .flip_plot(plot_out, flipped, add_x_text)
    return(plot_out)
}

#' @importFrom ggplot2 ggplot aes labs geom_point geom_raster
.feature_plotter <- function(
        feature_data,
        name,
        xlab = "Samples",
        flipped,
        add_legend,
        add_x_text,
        point_alpha,
        point_size){
    # If the values are factors, use coloring to plot them. This step is to
    # ensure that this functions works both with factors and numeric values.
    if(is.factor(feature_data$Y)){
        feature_data$colour_by <- feature_data$Y
        feature_data$Y <- ""
        colour_by <- unique(feature_data$feature_name)
    }
    # Start plotting
    feature_plot_out <- ggplot(
        feature_data, aes(x=.data[["X"]], y=.data[["Y"]])) +
        labs(x = xlab, y = name)
    # If there is only one value, i.e., the variable to be plotted was factor
    if(length(unique(feature_data$Y)) == 1L){
        # Create a bar layout
        feature_out <- .get_bar_args(
            colour_by, alpha = point_alpha, add_border = FALSE)
        feature_plot_out <- feature_plot_out +
            do.call(geom_raster, feature_out$args) +
            scale_y_discrete(expand = c(0,0))
        # Adjust the colour scale and legend title
        feature_plot_out <- .resolve_plot_colours(
            feature_plot_out, feature_data$colour_by, colour_by, fill = TRUE)
        legend_pos <- "bottom"
    } else {
        # If the values are numeric, create a point layout
        feature_out <- .get_point_args(
            NULL, shape_by = NULL, size_by = NULL, alpha = point_alpha,
            size = point_size)
        feature_plot_out <- feature_plot_out +
            do.call(geom_point, feature_out$args)
        legend_pos <- "right"
    }
    # Adjust theme
    feature_plot_out <- feature_plot_out +
        theme_classic()
    # Remove legend if specified, adjust the position
    feature_plot_out <- .add_legend(feature_plot_out, add_legend, legend_pos)
    # Flip the plot if specified
    feature_plot_out <- .flip_plot(feature_plot_out, flipped, add_x_text)
    return(feature_plot_out)
}

.features_plotter <- function(
        features_data,
        order_sample_by,
        xlab = NULL,
        flipped = FALSE,
        add_legend = add.legend,
        add.legend = TRUE,
        add_x_text = add.x.text,
        add.x.text = FALSE,
        point_alpha = point.alpha,
        point.alpha = 1,
        point_size = point.size,
        point.size = 2,
        ...){
    # Get the name of sample metadata variables that will be plotted
    names <- colnames(features_data)
    # For each variable, create a data.frame that contains sample names,
    # variable name and values of variable
    features_data <- lapply(
        names, function(col){
            data.frame(
                X = factor(rownames(features_data), rownames(features_data)),
                feature_name = col,
                Y = features_data[[col]])
        })
    names(features_data) <- names
    # Loop through variables and create plot for each variable
    plots_out <- mapply(
        .feature_plotter,
        features_data,
        names(features_data),
        MoreArgs = list(
            xlab = xlab,
            flipped = flipped,
            add_legend = add_legend,
            add_x_text = add_x_text,
            point_alpha = point_alpha,
            point_size = point_size),
        SIMPLIFY = FALSE)
    names(plots_out) <-  names(features_data)
    # If the varoable for order the data was specified, return only the feature
    # plot with that variable along with the main plot. This means that all
    # other feature plots are discarded.
    if(!is.null(order_sample_by)){
        reorder <- c(
            order_sample_by,
            names(plots_out)[!(names(plots_out) %in% order_sample_by)])
        m <- match(reorder,names(plots_out))
        m <- m[!is.na(m)]
        plots_out <- plots_out[m]
    }
    return(plots_out)
}
