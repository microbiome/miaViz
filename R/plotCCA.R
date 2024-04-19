#' Plot RDA or CCA object
#'
#' \code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
#' output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
#' for supervised ordination of microbiome data.
#'
#' @param object a
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#'   or a matrix of weights. The latter is returned as output from \code{\link[mia:runCCA]{getRDA}}.
#' 
#' @param dimred A string or integer scalar indicating the reduced dimension to
#'   plot. This is the output of \code{\link[mia:runCCA]{runRDA}} and resides in
#'   \code{reducedDim(tse, dimred)}.
#' 
#' @param add.ellipse One of \code{c(TRUE, FALSE, "fill", "colour", "color")},
#'   indicating whether ellipses should be present, absent, filled or colored.
#'   (default: \code{ellipse.fill = TRUE})
#'
#' @param ellipse.alpha Number between 0 and 1 to adjust the opacity of ellipses.
#'   (default: \code{ellipse.alpha = 0.2})
#'
#' @param ellipse.linewidth Number specifying the size of ellipses.
#'   (default: \code{ellipse.linewidth = 0.1})
#' 
#' @param ellipse.linetype Discrete number specifying the style of ellipses.
#'   (default: \code{ellipse.linetype = 1})
#'
#' @param add.vectors TRUE or FALSE, should vectors appear in the plot.
#'    (default: \code{add.vectors = TRUE})
#' 
#' @param vec.size Number specifying the size of vectors.
#'   (default: \code{vec.size = 0.5})
#' 
#' @param vec.colour String specifying the colour of vectors.
#'   (default: \code{vec.color = "black"})
#' 
#' @param vec.color Alias for `vec.colour`.
#' 
#' @param vec.linetype Discrete number specifying the style of vector lines.
#'   (default: \code{vec.linetype = 1})
#' 
#' @param arrow.size Number specifying the size of arrows.
#'   (defaults: \code{arrow.size = 0.25})
#' 
#' @param label.size Number specifying the size of text and labels.
#'   (default: \code{label.size = 4})
#' 
#' @param label.colour String specifying the colour of text and labels.
#'   (default: \code{label.color = "black"})
#' 
#' @param label.color Alias for `label.colour`.
#' 
#' @param sep.group String specifying the separator used in the labels.
#'   (default: \code{sep.group = "\U2014"})
#'   
#' @param repl.underscore String used to replace underscores in the labels.
#'   (default: \code{repl.underscore = " "})
#' 
#' @param vec.text TRUE or FALSE, should text instead of labels be used to label vectors.
#'   (default: \code{vec.text = TRUE})
#' 
#' @param repel.labels TRUE or FALSE, should labels be repelled.
#'   (default: \code{repel.labels = TRUE})
#'
#' @param parse.labels TRUE or FALSE, should labels be parsed.
#'   (default: \code{parse.labels = TRUE})
#'
#' @param add.significance TRUE or FALSE, should explained variance and p-value
#'   appear in the labels. (default: \code{add.significance = TRUE})
#'
#' @param add.expl.var TRUE or FALSE, should explained variance appear on the
#'   coordinate axes. (default: \code{add.expl.var = TRUE})
#' 
#' @param ... additional parameters for plotting, inherited from
#'   \code{\link[scater:plotReducedDim]{plotReducedDim}},
#'   \code{\link[ggplot2:geom_label]{geom_label}} and
#'   \code{\link[ggrepel:geom_label_repel]{geom_label_repel}}.
#' 
#' @details
#' \code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
#' output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
#' for supervised ordination of microbiome data. Either a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' or a matrix object is supported as input. When the input is a
#' TreeSummarizedExperiment, this should contain the output of runRDA
#' in the reducedDim slot and the argument \code{dimred} needs to be defined.
#' When the input is a matrix, this should be returned as output from
#' getRDA. However, the first method is recommended because it provides
#' the option to adjust aesthetics to the colData variables through the
#' arguments inherited from \code{\link[scater:plotReducedDim]{plotReducedDim}}.
#' 
#' @return 
#' A \code{ggplot2} object 
#'
#' @name plotCCA
#'
#' @examples
#' # Load dataset
#' library(miaViz)
#' data("enterotype", package = "mia")
#' tse <- enterotype
#'  
#' # Run RDA and store results into TreeSE
#' tse <- runRDA(tse,
#'               formula = assay ~ ClinicalStatus + Gender + Age,
#'               FUN = vegan::vegdist,
#'               distance = "bray",
#'               na.action = na.exclude)
#'                
#' # Create RDA plot coloured by variable
#' plotRDA(tse, "RDA",
#'         colour_by = "ClinicalStatus")
#'  
#' # Create RDA plot with empty ellipses
#' plotRDA(tse, "RDA",
#'         colour_by = "ClinicalStatus",
#'         add.ellipse = "colour")
#'  
#' # Create RDA plot with text encased in labels
#' plotRDA(tse, "RDA",
#'         colour_by = "ClinicalStatus",
#'         vec.text = FALSE)
#'  
#' # Create RDA plot without repelling text
#' plotRDA(tse, "RDA",
#'         colour_by = "ClinicalStatus",
#'         repel.labels = FALSE)
#'  
#' # Create RDA plot without vectors
#' plotRDA(tse, "RDA",
#'         colour_by = "ClinicalStatus",
#'         add.vectors = FALSE)
#'  
#' # Calculate RDA as a separate object
#' rda_mat <- getRDA(tse,
#'                         formula = assay ~ ClinicalStatus + Gender + Age,
#'                         FUN = vegan::vegdist,
#'                         distance = "bray",
#'                         na.action = na.exclude)
#'  
#' # Create RDA plot from RDA matrix
#' plotRDA(rda_mat)
NULL

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setGeneric("plotCCA", signature = c("object"),
           function(object, ...) standardGeneric("plotCCA"))

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setMethod("plotCCA", signature = c(object = "SingleCellExperiment"),
    function(object, dimred, ...){
        # Reproduce plotRDA function
        return(plotRDA(object, dimred, ...))
    }
)

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setMethod("plotCCA", signature = c(object = "matrix"),
    function(object, ...){
        # Reproduce plotRDA function
        return(plotRDA(object, ...))
    }
)

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setGeneric("plotRDA", signature = c("object"),
    function(object, ...) standardGeneric("plotRDA"))

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setMethod("plotRDA", signature = c(object = "SingleCellExperiment"),
    function(object, dimred,
             add.ellipse = TRUE, ellipse.alpha = 0.2, ellipse.linewidth = 0.1, ellipse.linetype = 1,
             vec.size = 0.5, vec.color = vec.colour, vec.colour = "black", vec.linetype = 1,
             arrow.size = 0.25, label.color = label.colour, label.colour = "black", label.size = 4,
             vec.text = TRUE, repel.labels = TRUE, sep.group = "\U2014", repl.underscore = " ",
             add.significance = TRUE, add.expl.var = TRUE, add.vectors = TRUE, parse.labels = TRUE, ...){
        ###################### Input check ########################
        if( !(add.ellipse %in% c(TRUE, FALSE, "fill", "color", "colour")) ){
            stop("'add.ellipse' must be one of c(TRUE, FALSE, 'fill', 'color', 'colour').", call. = FALSE)
        }
        if ( !.is_a_bool(add.vectors) ){
            stop("'add.vectors must be TRUE or FALSE.", call. = FALSE)
        }
        if ( !add.vectors ){
            warning("'add.vectors' is FALSE, so other arguments for vectors and labels will be disregarded.", call. = FALSE)
        }
        if( !.is_a_bool(vec.text) ){
            stop("'vec.text' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(repel.labels) ){
            stop("'repel.labels' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(parse.labels) ){
            stop("'parse.labels' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(add.significance) ){
            stop("'add.significance' must be TRUE or FALSE.", call. = FALSE)
        }
        if( parse.labels && !add.significance ){
            parse.labels <- FALSE
            warning("'parse.labels' was turned off because 'add.significance' is FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(add.expl.var) ){
            stop("'add.expl.var' must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.are_whole_numbers(ellipse.linetype) ){
            stop("'ellipse.linetype' must be a whole number.", call. = FALSE)
        }
        if( !.are_whole_numbers(ellipse.linetype) ){
            stop("'vec.linetype' must be a whole number.", call. = FALSE)
        }
        if( ellipse.alpha < 0 || ellipse.alpha > 1 ){
            stop("'ellipse.alpha' must be a number between 0 and 1.", call. = FALSE)
        }
        if ( !is.numeric(ellipse.linewidth) && ellipse.linewidth > 0 ) {
            stop("'ellipse.linewidth' must be a positive number.", call. = FALSE)
        }
        if ( !is.numeric(vec.size) && vec.size > 0 ) {
            stop("'vec.size' must be a positive number.", call. = FALSE)
        }
        if ( !is.numeric(arrow.size) && arrow.size > 0 ) {
            stop("'arrow.size' must be a positive number.", call. = FALSE)
        }
        if ( !is.numeric(label.size) && label.size > 0 ) {
            stop("'label.size' must be a positive number.", call. = FALSE)
        }
        if ( !.is_non_empty_string(vec.color) ) {
            stop("'vec.color' must be a non-empty string specifying a colour", call. = FALSE)
        }
        if ( !.is_non_empty_string(label.color) ) {
            stop("'label.color' must be a non-empty string specifying a colour", call. = FALSE)
        }
        if ( !.is_a_string(sep.group) ) {
            stop("'sep.group' must be a string specifying a separator.", call. = FALSE)
        }
        if ( !.is_a_string(repl.underscore) ) {
            stop("'repl.underscore' must be a string.", call. = FALSE)
        }
        ###################### Input check end ####################
        # Get data for plotting
        plot_data <- .incorporate_rda_vis(
            object, dimred, sep.group = sep.group, repl.underscore = repl.underscore,
            add.significance = add.significance, add.expl.var = add.expl.var,
            add.ellipse = add.ellipse, add.vectors = add.vectors, ...
        )
        # Create a plot
        plot <- .rda_plotter(
            plot_data, ellipse.alpha = ellipse.alpha, ellipse.linewidth = ellipse.linewidth,
            ellipse.linetype = ellipse.linetype, vec.size = vec.size, vec.color = vec.color,
            vec.colour = vec.colour, vec.linetype = vec.linetype, arrow.size = arrow.size,
            label.color = label.color, label.colour = label.colour, label.size = label.size,
            vec.text = vec.text, add.ellipse = add.ellipse, repel.labels = repel.labels,
            parse.labels = parse.labels, ...
        )
        return(plot)
    }
)

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setMethod("plotRDA", signature = c(object = "matrix"),
    function(object, ...){
        # Construct TreeSE from rda/cca object
        object <- .rda2tse(object)
        # Run plotRDA method for TreeSE
        return(plotRDA(object, "RDA", ...))
    }
)


################## HELP FUNCTIONS ##########################

# Construct TreeSE from rda/cca object to pass it to downstream functions
.rda2tse <- function(object) {
    # Convert rda/cca object to TreeSE
    object <- TreeSummarizedExperiment(
      assays = matrix(ncol = nrow(object), dimnames = list(NULL, rownames(object))),
      reducedDims = list(RDA = object)
    )
    return(object)
}

# Get data for plotting
#' @importFrom scater plotReducedDim retrieveCellInfo
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
.incorporate_rda_vis <- function(
        tse, dimred, ncomponents = 2, colour_by = color_by, color_by = NULL,
        add.significance = TRUE, add.expl.var = TRUE, add.ellipse = TRUE,
        add.vectors = TRUE, vec.lab = NULL, expl_var = NULL,
        sep.group = "\U2012", repl.underscore = " ", ...){

    # Check dimred
    if( !dimred %in% reducedDimNames(tse) ){
        stop("'dimred' must specify reducedDim.", call. = FALSE)
    }
    # Get reducedDim
    reduced_dim <- reducedDim(tse, dimred)
    # Check that there are at least 2 coordinates
    if( ncol(reduced_dim) < 2 ){
        stop("reducedDim specified by 'dimred' must have at least 2 columns.", call. = FALSE)
    }
    # Only 2 dimensions are supported currently
    ncomponents <- 2
    
    # If specified, get explained variance
    if( add.expl.var ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("rda", "cca")
        # If it can be found
        if( any(ind) ){
            # Add explained variance
            rda <- attributes(reduced_dim)[ind][[1]]
            expl_var <- summary(rda)$concont$importance[2, ]*100
        } else{
            # If it cannot be found, give warning
            warning(paste("RDA/CCA object was not found. Please compute",
                          "RDA/CCA by using runCCA or getCCA."),
                    call. = FALSE)
        }
    }
    
    # Get scatter plot with plotReducedDim --> keep theme similar between ordination methods
    plot <- plotReducedDim(
        tse, dimred = dimred, ncomponents = ncomponents, colour_by = colour_by,
        percentVar = expl_var, ...)
    
    # Get data for ellipse
    ellipse_data <- NULL
    if( add.ellipse != FALSE && !is.null(colour_by) ){
        ellipse_data <- reduced_dim
        ellipse_data <- as.data.frame(ellipse_data)
        ellipse_data[[colour_by]] <- retrieveCellInfo(tse, colour_by)[["value"]]
        attributes(ellipse_data)$colour_by <- colour_by
    }
    
    # Get data for vectors
    vector_data <- NULL
    if( add.vectors ){
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
            warning(paste("RDA/CCA object was not found. Please compute RDA/CCA",
                          "by using runCCA or getCCA."),
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
        if( is.null(vec.lab) ){
            vector_label <- rownames(vector_data)
            # Make labels more tidy
            if( all_var_found ){
                vector_label <- .tidy_vector_labels(
                    vector_label, coldata,
                    sep.group = sep.group,
                    repl.underscore = repl.underscore
                )
            }
            # Add to df
            vector_data$vector_label <- vector_label
        } else{
            # Check that user-provided labels are correct length
            if( length(vec.lab) != nrow(vector_data) ){
                stop("Number of labels in 'vec_lab' do not match with number of vectors.",
                     call. = FALSE)
            }
            # If they are, add labels to data
            vector_data$vector_label <- vec.lab
        }
    }
    
    # Get significance data
    signif_data <- NULL
    if( add.significance && !is.null(vector_data) && all_var_found ){
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
            vector_label <- .add_signif_to_vector_labels(
                vector_label, variable_names, signif_data, repl.underscore)
            vector_data[["vector_label"]] <- vector_label
        } else{
            # If it cannot be found, give warning
            warning(paste("Significance data was not found. please compute",
                          "CCA/RDA by using runCCA or getCCA."),
                    call. = FALSE)
        }
    }
    
    # Create a list to return
    result <- list(
        plot = plot,
        ellipse_data = ellipse_data,
        vector_data = vector_data
    )
    return(result)
}

# Make vector labels more tidy, i.e, separate variable and group names.
# Replace also underscores with space
.tidy_vector_labels <- function(
        vector_label, coldata, sep.group = "\U2014", repl.underscore = " ", ...){
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
            new_name <- paste0(var_name, " ", sep.group, " ", group_name)
        } else{
            new_name <- name
        }
        # Replace underscores with space
        new_name <- gsub("_", repl.underscore, new_name)
        return(new_name)
    })
    return(vector_label)
}

# This function adds significance info to vector labels
.add_signif_to_vector_labels <- function(
        vector_label, var_names, signif_data, repl.underscore = " ", ...){
    # Replace underscores from significance data and variable names to match labels
    rownames(signif_data) <- sapply(rownames(signif_data), function(x) gsub("_", repl.underscore, x))
    var_names <- sapply(var_names, function(x) gsub("_", repl.underscore, x))
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
        plot_data, ellipse.alpha = 0.2, ellipse.linewidth = 0.1, ellipse.linetype = 1,
        vec.size = 0.5, vec.color = vec.colour, vec.colour = "black",
        vec.linetype = 1, arrow.size = 0.25, min.segment.length = 5,
        label.color = label.colour, label.colour = "black", label.size = 4,
        parse.labels = TRUE, vec.text = TRUE, repel.labels = TRUE, add.ellipse = TRUE,
        position = NULL, nudge_x = NULL, nudge_y = NULL, direction = "both",
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
        if( add.ellipse %in% c(TRUE, "fill") ){
            plot <- plot +
                stat_ellipse(data = data,
                             aes(x = .data[[xvar]], y = .data[[yvar]],
                                 color = .data[[colour_var]], fill = after_scale(color)),
                             geom = "polygon", alpha = ellipse.alpha,
                             linewidth = ellipse.linewidth, linetype = ellipse.linetype)
        } else if ( add.ellipse %in% c("color", "colour") ){
            plot <- plot +
                stat_ellipse(data = data,
                             aes(x = .data[[xvar]], y = .data[[yvar]], color = .data[[colour_var]]),
                             geom = "polygon", alpha = 0, linetype = ellipse.linetype)
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
                         arrow = arrow(length = unit(arrow.size, "cm")),
                         color = vec.color, linetype = vec.linetype, size = vec.size)
        # Add vector labels (text or label)
        # Make list of arguments for geom_text/geom_label
        label_args <- list(
                        data = data,
                        mapping = aes(x = .data[[xvar]], y = .data[[yvar]]),
                        label = data[["vector_label"]], parse = parse.labels,
                        color = label.color, size = label.size, stat = "identity",
                        nudge_x = nudge_x, nudge_y = nudge_y, show.legend = NA,
                        na.rm = FALSE, inherit.aes = TRUE
                      )
        # Repel text
        if( repel.labels ){
          # Add arguments for geom_text_repel/geom_label_repel to list
          label_args <- c(
            label_args, min.segment.length = min.segment.length,
            box.padding = 0.25, point.padding = 1e-06, force = 1, force_pull = 1,
            max.time = 0.5, max.iter = 10000, max.overlaps = max.overlaps,
            direction = direction, seed = NA, verbose = FALSE
          )
          # repelled text
          if( vec.text ){
            plot <- plot + do.call("geom_text_repel", label_args)
          # repelled labels
          } else{
            plot <- plot + do.call("geom_label_repel", label_args)
          }
        # Do not repel text
        } else{
          # not repelled text
          if( vec.text ){
            label_args <- c(label_args, check_overlap = check_overlap)
            plot <- plot + do.call("geom_text", label_args)
          # not repelled labels
          } else{
            plot <- plot + do.call("geom_label", label_args)
          }
        }
        
    }
    return(plot)
}
