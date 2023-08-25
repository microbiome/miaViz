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
#' \dontrun{
#' library(mia)
#' 
#' tse <- runRDA()
#' plotRDA()
#' 
#' 
#' rda <- calculateCCA()
#' plotCCA()
#' 
#' rda <- attributes(tse)$rda
#' 
#' plotRDA()
NULL

#' @rdname plotRDA
#' @export
setGeneric("plotRDA", signature = c("object"),
    function(object, ...) standardGeneric("plotRDA"))


#' @rdname plotRDA
#' @export
setMethod("plotRDA", signature = c(object = "SingleCellExperiment"),
    function(object, ...){
        ###################### Input check #######################
        
        ###################### Input check end ####################
        # Get data for plotting
        plot_data <- .incorporate_rda_vis(object, ...)
        # Create a plot
        plot <- .rda_plotter(plot_data, ...)
        return(plot)
    }
)

################## HELP FUNCTIONS ##########################
# Get data for plotting
.incorporate_rda_vis <- function(
        tse, dimred, ncomponents = 2, colour_by = color_by, color_by = NULL,
        add_ellipse = TRUE, add_vectors = TRUE, add_significance = TRUE, add_expl_var = TRUE,
        vec_lab = NULL, ...){
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
    
    # Get scatter plot with plotReducedDim --> keep theme similar between ordination methods
    plot <- plotReducedDim(tse, dimred = dimred, ncomponents = ncomponents, colour_by = colour_by, ...)
    
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
            warning("CCA object was not found. Please calculate CCA by using runCCA/calculateCCA.", call. = FALSE)
        }
    }
    # Get variable names from sample metadata
    variable_names <- colnames(colData(tse))
    # Get vector labels
    if( !is.null(vector_data) ){
        # If user has not provided vector labels
        if( is.null(vec_lab) ){
            vector_label <- rownames(vector_data)
            # Make labels more tidy
            if( length(variable_names) > 0 ){
                vector_label <- .tidy_vector_labels(vector_label, variable_names)
            }
            # Add to df
            vector_data$vector_label <- vector_label
        } else{
            # Check that user-provided labels are correct length
            if( length(vec_lab) != nrow(vector_data) ){
                stop("Number of labels in 'vec_lab' do not match with number of vectors.", call. = FALSE)
            }
            # If they are, add labels to data
            vector_data$vector_label <- vec_lab
        }
    }
    
    # Get significance data
    signif_data <- NULL
    if( add_significance && !is.null(vector_data) ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("significance")
        # If it can be found
        if( any(ind) ){
            # Get biplot from cca object
            signif_data <- attributes(reduced_dim)[ind][[1]]
            signif_data <- signif_data[["permanova"]]
            signif_data <- as.data.frame(signif_data)
        } else{
            # If it cannot be found, give warning
            warning("Significance data was not found. please calculate CCA by using runCCA/calculateCCA.", call. = FALSE)
        }
    }
    # Add significance to vector labels
    if( !is.null(signif_data) ){
        # Get vector labels
        vector_label <- vector_data[["vector_label"]]
        # Add significance to vector labels
        if( length(variable_names) > 0 ){
            vector_label <- .add_signif_to_vector_labels(vector_label, variable_names)
        }
        vector_data[["vector_label"]] <- vector_label
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
            warning("CCA object was not found. Please calculate CCA by using runCCA/calculateCCA.", call. = FALSE)
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

# Make vector labels more tidy
.tidy_vector_labels <- function(vector_label, variable_names){
    # Check that every label can be found from the variable names
    all_found <- all(colSums(
        sapply(rownames(rda$CCA$biplot), function(x)
            sapply(variable_names, function(y) grepl(y, x)) )) == 1)
    if( all_found ){
        # Loop through vector labels
        vector_label <- sapply(vector_label, FUN = function(name){
            # Get the real variable name from sample metadata
            variable_name <- variable_names[ sapply(variable_names, function(x) grepl(x, name)) ]
            # If the vector label includes also group name
            if( !name %in% variable_names ){
                # Get the group name
                group_name <- unique( coldata[[variable_name]] )[ 
                    which( paste0(variable_name, unique( coldata[[variable_name]] )) == name ) ]
                # Modify vector so that group is separated from variable name
                new_name <- paste0(variable_name, " \U2012 ", group_name)
            } else{
                new_name <- name
            }
            return(new_name)
        })
    }
    return(vector_label)
}

.add_signif_to_vector_labels <- function(vector_label, variable_names){
    # Check that every label can be found from the variable names
    all_found <- all(colSums(
        sapply(rownames(rda$CCA$biplot), function(x)
            sapply(variable_names, function(y) grepl(y, x)) )) == 1)
    if( all_found ){
        # Loop through vector labels
        vector_label <- sapply(vector_label, FUN = function(name){
            # Get the real variable name from sample metadata
            variable_name <- variable_names[ sapply(variable_names, function(x) grepl(x, name)) ]
            # Add percentage how much this variable explains, and p-value
            new_name <- expr(
                paste(!!name, " (", 
                      !!format(
                          round( signif_data[variable_name, "Explained variance"]*100, 1),
                          nsmall = 1), "%, ", italic("P"), " = ",
                      !!gsub("0\\.","\\.", format(
                          round( signif_data[variable_name, "Pr(>F)"], 3),
                          nsmall = 3)), ")"))
            
            return(new_name)
        })
    }
    return(vector_label)
}

# Plot based on the data
.rda_plotter <- function(plot_data, alpha = 0.2, vec_size = 0.25, vec_color = vec_colour, vec_colour = "black", ...){
    # Get the scatter plot
    plot <- plot_data[["plot"]]
    # Add ellipse
    if( !is.null(plot_data$ellipse_data) ){
        # Get data and variabe names
        data <- plot_data$ellipse_data
        xvar <- colnames(data)[[1]]
        yvar <- colnames(data)[[2]]
        colour_var <- attributes(plot_data$ellipse_data)$colour_by
        # Add ellipses to plot
        plot <- plot +
            stat_ellipse(data = data,
                         aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[colour_var]]),
                         geom="polygon", alpha = alpha)
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
                         arrow = arrow(length = unit(vec_size, "cm")), color = vec_color)
        # Add vector labels
        plot <- plot +
            geom_text_repel(data = data, aes(x = .data[[xvar]], y = .data[[yvar]]),
                            label = data[["vector_label"]], parse = TRUE
            )
        
    }
    # Add axis labels
    plot <- plot + xlab(plot_data$xlab) + ylab(plot_data$ylab)
    return(plot)
}
