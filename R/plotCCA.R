#' Plot RDA or CCA object
#'
#'
#' @param object a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'    object or (db)RDA or CCA object..
#'
#' @param assay.type a single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} to be
#'   plotted. (default: \code{assay.type = "counts"})
#'   
#'   
#' @param ... additional parameters for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}
#'
#' @details
#' This function creates a RDA/CCA plot. If the input is TreeSE, the input
#' 
#' If the input is result of calculateRDA...
#' 
#' If the input is rda/cca object...
#' 
#' 
#' @return 
#' A \code{ggplot2} object 
#'
#' @name plotRDA
#'
#' @examples
#' 
#'  data(peerj13075)
#'  tse <- peerj13075
#'  # Perform RDA
#'  tse <- runRDA(tse, data ~ Age + Geographical_location, FUN = vegan::vegdist, method = "bray")
#'  # Create a plot
#'  plotRDA(tse, "RDA", colour_by = "Geographical_location")
#'  
NULL

#' @rdname plotRDA
#' @export
setGeneric("plotRDA", signature = c("object"),
    function(object, dimred, ...) standardGeneric("plotRDA"))


#' @rdname plotRDA
#' @export
setMethod("plotRDA", signature = c(object = "SingleCellExperiment"),
    function(object, dimred, ...){
        ###################### Input check #######################
        
        ###################### Input check end ####################
        # Get data for plotting
        plot_data <- .incorporate_rda_vis(object, dimred, ...)
        # Create a plot
        plot <- .rda_plotter(plot_data, ...)
        return(plot)
    }
)

################## HELP FUNCTIONS ##########################

# Get data for plotting
#' @importFrom scater plotReducedDim
.incorporate_rda_vis <- function(
        tse, dimred, ncomponents = 2, colour_by = color_by, color_by = NULL,
        shape_by = NULL, size_by = NULL, order_by = NULL, text_by = NULL,
        other_fields = list(), swap_rownames = NULL, point.padding = NA,
        add_ellipse = TRUE, add_vectors = TRUE, add_significance = TRUE,
        add_expl_var = TRUE, vec_lab = NULL, bins = NULL, ...){
    
    # Check dimred
    if( !dimred %in% reducedDimNames(tse) ){
        stop("'dimred' must specify reducedDim.", call. = FALSE)
    }
    # Get reducedDim
    reduced_dim <- reducedDim(tse, dimred)
    # Check that there are at least 2 coordinates.
    if( ncol(reduced_dim) < 2 ){
        stop("reducedDim specified by 'dimred' must have at least 2 columns.", call. = FALSE)
    }
    # Only 2 dimensions are supported currently
    ncomponents <- 2
    # Make list of arguments for plotReducedDim
    plotReducedDim_args <- list(
      object = tse, dimred = dimred, ncomponents = ncomponents, colour_by = color_by,
      color_by = color_by, shape_by = shape_by, size_by = size_by, order_by = order_by,
      text_by = text_by, text_size = 5, text_colour = "black", text_color = "black",
      label_format = c("%s %i", " (%i%%)"), other_fields = other_fields,
      swap_rownames = swap_rownames, point.padding = point.padding, force = 1,
      rasterise = FALSE, scattermore = FALSE, summary_fun = "sum", hex = FALSE
    )
    # Get scatter plot with plotReducedDim --> keep theme similar between ordination methods
    plot <- do.call("plotReducedDim", plotReducedDim_args)
    
    # Get data for ellipse
    ellipse_data <- NULL
    if( add_ellipse && !is.null(colour_by) ){
        ellipse_data <- reduced_dim
        ellipse_data <- as.data.frame(ellipse_data)
        ellipse_data[[colour_by]] <- colData(tse)[[colour_by]]
        attributes(ellipse_data)$colour_by <- colour_by
    }
    
    # Get data for vectors
    vector_data <- NULL
    if( add_vectors ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("rda", "cca")
        # If it can be found
        if( any(ind) ){
            # Get biplot from cca object
            rda <- attributes(reduced_dim)[ind][[1]]
            vector_data <- rda$CCA$biplot
            vector_data <- as.data.frame(vector_data)
            vector_data[["group"]] <- rownames(vector_data)
        } else{
            # If it cannot be found, give warning
            warning(
                "CCA object was not found. Please calculate CCA by using runCCA.",
                call. = FALSE)
        }
    }
    
    # Get variable names from sample metadata
    coldata <- colData(tse)
    variable_names <- colnames(coldata)
    # Check if variable names can be found metadata
    all_var_found <- FALSE
    if( !is.null(vector_data) && length(variable_names) > 0 ){
        all_var_found <- all(colSums(
            sapply(rownames(vector_data), function(x)
                sapply(variable_names, function(y) grepl(y, x)) )) == 1)
    }
    
    # Get vector labels
    if( !is.null(vector_data) ){
        # Get those names that are present in data
        
        # If user has not provided vector labels
        if( is.null(vec_lab) ){
            vector_label <- rownames(vector_data)
            # Make labels more tidy
            if( all_var_found ){
                vector_label <- .tidy_vector_labels(vector_label, coldata, ...)
            }
            # Add to df
            vector_data$vector_label <- vector_label
        } else{
            # Check that user-provided labels are correct length
            if( length(vec_lab) != nrow(vector_data) ){
                stop("Number of labels in 'vec_lab' do not match with number of vectors.",
                     call. = FALSE)
            }
            # If they are, add labels to data
            vector_data$vector_label <- vec_lab
        }
    }
    
    # Get significance data
    signif_data <- NULL
    if( add_significance && !is.null(vector_data) && all_var_found ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("significance")
        # If it can be found
        if( any(ind) ){
            # Get biplot from cca object
            signif_data <- attributes(reduced_dim)[ind][[1]]
            signif_data <- signif_data[["permanova"]]
            signif_data <- as.data.frame(signif_data)
            # Add significance to vector labels
            # Get vector labels
            vector_label <- vector_data[["vector_label"]]
            # Add significance to vector labels
            vector_label <- .add_signif_to_vector_labels(vector_label, variable_names, signif_data, ...)
            vector_data[["vector_label"]] <- vector_label
        } else{
            # If it cannot be found, give warning
            warning("Significance data was not found. please calculate CCA by using runCCA.",
                    call. = FALSE)
        }
    }
    
    # Create labels for axis
    xlab <- paste0(dimred, " 1")
    ylab <- paste0(dimred, " 2")
    if( add_expl_var ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("rda", "cca")
        # If it can be found
        if( any(ind) ){
            # Add explained variance
            rda <- attributes(reduced_dim)[ind][[1]]
            xlab <- paste0(
                xlab, " (",
                format(round( summary(rda)$concont$importance[2, 1]*100, 1 ), nsmall = 1 ), "%)")
            ylab <- paste0(
                ylab, " (",
                format(round( summary(rda)$concont$importance[2, 2]*100, 1 ), nsmall = 1 ), "%)")
        } else{
            # If it cannot be found, give warning
            warning("CCA object was not found. Please calculate CCA by using runCCA.", call. = FALSE)
        }
    }
    
    # Create a list to return
    result <- list(
        plot = plot,
        ellipse_data = ellipse_data,
        vector_data = vector_data,
        xlab = xlab,
        ylab = ylab
    )
    return(result)
}

# Make vector labels more tidy, i.e, separate variable and group names.
# Replace also underscores with space
.tidy_vector_labels <- function(
        vector_label, coldata, sep_group = "\U2012", sep_underscore = " ", ...){
    # Get variable names from sample metadata
    var_names <- colnames(coldata)
    # Loop through vector labels
    vector_label <- sapply(vector_label, FUN = function(name){
        # Get the real variable name from sample metadata
        var_name <- var_names[ sapply(var_names, function(x) grepl(x, name)) ]
        # If the vector label includes also group name
        if( !name %in% var_names ){
            # Get the group name
            group_name <- unique( coldata[[var_name]] )[ 
                which( paste0(var_name, unique( coldata[[var_name]] )) == name ) ]
            # Modify vector so that group is separated from variable name
            new_name <- paste0(var_name, " ", sep_group, " ", group_name)
        } else{
            new_name <- name
        }
        # Replace underscores with space
        new_name <- gsub("_", sep_underscore, new_name)
        return(new_name)
    })
    return(vector_label)
}

# This function adds significance info to vector labels
.add_signif_to_vector_labels <- function(vector_label, var_names, signif_data, sep_underscore = " ", ...){
    # Replace underscores from significance data and variable names to match labels
    rownames(signif_data) <- sapply(rownames(signif_data), function(x) gsub("_", sep_underscore, x))
    var_names <- sapply(var_names, function(x) gsub("_", sep_underscore, x))
    # Loop through vector labels
    vector_label <- sapply(vector_label, FUN = function(name){
        # Get the real variable name from sample metadata
        var_name <- var_names[ sapply(var_names, function(x) grepl(x, name)) ]
        # Add percentage how much this variable explains, and p-value
        new_name <- expr(
            paste(!!name, " (", 
                  !!format(
                      round( signif_data[var_name, "Explained variance"]*100, 1),
                      nsmall = 1), "%, ", italic("P"), " = ",
                  !!gsub("0\\.","\\.", format(
                      round( signif_data[var_name, "Pr(>F)"], 3),
                      nsmall = 3)), ")"))
        
        return(new_name)
    })
    return(vector_label)
}

# Plot based on the data
#' @importFrom ggrepel geom_text_repel geom_label_repel
.rda_plotter <- function(
        plot_data, alpha = 0.2, ellipse_size = 0.1, ellipse_linetype = 1,
        vec_size = 0.5, vec_color = vec_colour, vec_colour = "black",
        vec_linetype = 1, arrow_size = 0.25, min.segment.length = 5,
        text_color = text_colour, text_colour = "black", text_size = 4,
        parse = TRUE, vec_text = TRUE, ellipse_fill = TRUE, repel_text = TRUE,
        nudge_x = 0, nudge_y = 0, direction = "both",
        max.overlaps = 10, check_overlap = FALSE, ...){

    # Get the scatter plot
    plot <- plot_data[["plot"]]
    # Add ellipse
    if( !is.null(plot_data$ellipse_data) ){
        # Get data and variabe names
        data <- plot_data$ellipse_data
        xvar <- colnames(data)[[1]]
        yvar <- colnames(data)[[2]]
        colour_var <- attributes(plot_data$ellipse_data)$colour_by
        # Add ellipses to plot (fill or colour the edge)
        if( ellipse_fill ){
            plot <- plot +
                stat_ellipse(data = data,
                             aes(x = .data[[xvar]], y = .data[[yvar]],
                                 color = .data[[colour_var]], fill = after_scale(color)),
                             geom = "polygon", alpha = alpha, size = ellipse_size,
                             linetype = ellipse_linetype)
        } else{
            plot <- plot +
                stat_ellipse(data = data,
                             aes(x = .data[[xvar]], y = .data[[yvar]], color = .data[[colour_var]]),
                             geom = "polygon", alpha = 0, linetype = ellipse_linetype)
        }
        
    }
    # Add vectors
    if( !is.null(plot_data$vector_data) ){
        # Get data and variabe names
        data <- plot_data$vector_data
        xvar <- colnames(data)[[1]]
        yvar <- colnames(data)[[2]]
        # Add vectors
        plot <- plot +
            geom_segment(data = data,
                         aes(x = 0, y = 0, xend = .data[[xvar]], yend = .data[[yvar]],
                             group = .data[["group"]]),
                         arrow = arrow(length = unit(arrow_size, "cm")),
                         color = vec_color, linetype = vec_linetype, size = vec_size)
        # Add vector labels (text or label)
        # Make list of arguments for geom_text/geom_label
        text_args <- list(
                      data = data,
                      mapping = aes(x = .data[[xvar]], y = .data[[yvar]]),
                      label = data[["vector_label"]], parse = parse,
                      color = text_color, size = text_size, stat = "identity",
                      nudge_x = nudge_x, nudge_y = nudge_y, show.legend = NA,
                      na.rm = FALSE, inherit.aes = TRUE
                     )
        # Repel text
        if( repel_text ){
          # Add arguments for geom_text_repel/geom_label_repel to list
          text_args <- c(
            text_args, min.segment.length = min.segment.length,
            box.padding = 0.25, point.padding = 1e-06, force = 1, force_pull = 1,
            max.time = 0.5, max.iter = 10000, max.overlaps = max.overlaps,
            direction = direction, seed = NA, verbose = FALSE
          )
          # repelled text
          if( vec_text ){
            plot <- plot + do.call("geom_text_repel", text_args)
          # repelled labels
          } else{
            plot <- plot + do.call("geom_label_repel", text_args)
          }
        # Do not repel text
        } else{
          # not repelled text
          if( vec_text ){
            text_args <- c(text_args, check_overlap = check_overlap)
            plot <- plot + do.call("geom_text", text_args)
          # not repelled labels
          } else{
            plot <- plot + do.call("geom_label", text_args)
          }
        }
        
    }
    # Add axis labels
    plot <- plot + xlab(plot_data$xlab) + ylab(plot_data$ylab)
    return(plot)
}
