#' Plot RDA or CCA object
#'
#' \code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
#' output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
#' for supervised ordination of microbiome data.
#'
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' or a matrix of weights. The latter is returned as output from
#' \code{\link[mia:runCCA]{getRDA}}.
#' 
#' @param dimred \code{Character scalar} or \code{integer scalar}. Determines
#' the reduced dimension to
#' plot. This is the output of \code{\link[mia:runCCA]{addRDA}} and resides in
#' \code{reducedDim(tse, dimred)}.
#' 
#' @param ... additional parameters for plotting, inherited from
#' \code{\link[scater:plotReducedDim]{plotReducedDim}},
#' \code{\link[ggplot2:geom_label]{geom_label}} and
#' \code{\link[ggrepel:geom_label_repel]{geom_label_repel}}.
#' \itemize{
#'   \item \code{add.ellipse}: One of
#'   \code{c(TRUE, FALSE, "fill", "colour")}, indicating whether
#'   ellipses should be present, absent, filled or colored.
#'   (default: \code{ellipse.fill = TRUE})
#'   
#'   \item \code{ellipse.alpha}: \code{Numeric scalar}. Between 0 and 1.
#'   Adjusts the opacity of ellipses. (Default: \code{0.2})
#'   
#'   \item \code{ellipse.linewidth}: \code{Numeric scalar}. Specifies the size
#'   of ellipses. (Default: \code{0.1})
#' 
#'   \item \code{ellipse.linetype}: \code{Integer scalar}. Specifies the style
#'   of ellipses. (Default: \code{1})
#' 
#'   \item \code{confidence.level}: \code{Numeric scalar}. Between 0 and 1.
#'   Adjusts confidence level. (Default: \code{0.95})
#' 
#'   \item \code{add.vectors}: \code{Logical scalar}. Should vectors appear in
#'   the plot. (Default: \code{TRUE})
#' 
#'   \item \code{vec.size}: \code{Numeric scalar}. Specifies the size of
#'   vectors. (Default: \code{0.5})
#' 
#'   \item \code{vec.colour}: \code{Character scalar}. Specifies the colour of
#'   vectors. (Default: \code{"black"})
#' 
#'   \item \code{vec.linetype}: \code{Integer scalar}. Specifies the style of
#'   vector lines. (Default: \code{1})
#' 
#'   \item \code{arrow.size}: \code{Numeric scalar}. Specifies the size of
#'   arrows. (Default: \code{arrow.size = 0.25})
#' 
#'   \item \code{label.size}: \code{Numeric scalar}. Specifies the size of text
#'   and labels. (Default: \code{4})
#' 
#'   \item \code{label.colour}: \code{Character scalar}. Specifies the colour of
#'   text and labels. (Default: \code{"black"})
#' 
#'   \item \code{sep.group}: \code{Character scalar}. Specifies the separator
#'   used in the labels. (Default: \code{"\U2014"})
#' 
#'   \item \code{repl.underscore}: \code{Character scalar}. Used to replace
#'   underscores in the labels. (Default: \code{" "})
#' 
#'   \item \code{vec.text}: \code{Logical scalar}. Should text instead of labels
#'   be used to label vectors. (Default: \code{TRUE})
#' 
#'   \item \code{repel.labels}: \code{Logical scalar}. Should labels be
#'   repelled. (Default: \code{TRUE})
#' 
#'   \item \code{parse.labels}: \code{Logical scalar}. Should labels be parsed.
#'   (Default: \code{TRUE})
#' 
#'   \item \code{add.significance}: \code{Logical scalar}. Should explained
#'   variance and p-value appear in the labels. (Default: \code{TRUE})
#' 
#'   \item \code{add.expl.var}: \code{Logical scalar}. Should explained
#'   variance appear on the coordinate axes. (Default: \code{TRUE})
#'   
#'   \item \code{add.centroids}: \code{Logical scalar}. Should centroids
#'   of variables be added. (Default: \code{FALSE})
#'   
#'   \item \code{add.species}: \code{Logical scalar}. Should species
#'   scores be added. (Default: \code{FALSE})
#' }
#'
#' 
#' @details
#' \code{plotRDA} and \code{plotCCA} create an RDA/CCA plot starting from the
#' output of \code{\link[mia:runCCA]{CCA and RDA}} functions, two common methods
#' for supervised ordination of microbiome data. Either a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-constructor]{TreeSummarizedExperiment}}
#' or a matrix object is supported as input. When the input is a
#' TreeSummarizedExperiment, this should contain the output of addRDA
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
#' tse <- addRDA(
#'     tse,
#'     formula = assay ~ ClinicalStatus + Gender + Age,
#'     FUN = getDissimilarity,
#'     distance = "bray",
#'     na.action = na.exclude
#'     )
#'                
#' # Create RDA plot coloured by variable
#' plotRDA(tse, "RDA", colour.by = "ClinicalStatus")
#'  
#' # Create RDA plot with empty ellipses
#' plotRDA(tse, "RDA", colour.by = "ClinicalStatus", add.ellipse = "colour")
#'  
#' # Create RDA plot with text encased in labels
#' plotRDA(tse, "RDA", colour.by = "ClinicalStatus", vec.text = FALSE)
#'  
#' # Create RDA plot without repelling text
#' plotRDA(tse, "RDA", colour.by = "ClinicalStatus", repel.labels = FALSE)
#'  
#' # Create RDA plot without vectors
#' plotRDA(tse, "RDA", colour.by = "ClinicalStatus", add.vectors = FALSE)
#'  
#' # Calculate RDA as a separate object
#' rda_mat <- getRDA(
#'     tse,
#'     formula = assay ~ ClinicalStatus + Gender + Age,
#'     FUN = getDissimilarity,
#'     distance = "bray",
#'     na.action = na.exclude
#'     )
#'  
#' # Create RDA plot from RDA matrix
#' plotRDA(rda_mat)
NULL

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setGeneric("plotCCA", signature = c("x"),
    function(x, ...) standardGeneric("plotCCA"))

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setMethod("plotCCA", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ...){
        # Reproduce plotRDA function
        return(plotRDA(x, dimred, ...))
    }
)

#' @rdname plotCCA
#' @aliases plotRDA
#' @export
setMethod("plotCCA", signature = c(x = "matrix"),
    function(x, ...){
        # Reproduce plotRDA function
        return(plotRDA(x, ...))
    }
)

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setGeneric("plotRDA", signature = c("x"),
    function(x, ...) standardGeneric("plotRDA"))

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setMethod("plotRDA", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ...){
        ###################### Input check ########################
        # Check dimred
        if( !dimred %in% reducedDimNames(x) ){
            stop("'dimred' must specify reducedDim.", call. = FALSE)
        }
        ###################### Input check end ####################
        # Get data for plotting
        plot_data <- .incorporate_rda_vis(x, dimred, ...)
        # Create a plot
        plot <- .rda_plotter(plot_data, ...)
        return(plot)
    }
)

#' @rdname plotCCA
#' @aliases plotCCA
#' @export
setMethod("plotRDA", signature = c(x = "matrix"),
    function(x, ...){
        # Construct TreeSE from rda/cca object
        x <- .rda2tse(x)
        # Run plotRDA method for TreeSE
        return(plotRDA(x, "RDA", ...))
    }
)


################## HELP FUNCTIONS ##########################

# Construct TreeSE from rda/cca object to pass it to downstream functions
.rda2tse <- function(object) {
    # Convert rda/cca object to TreeSE
    object <- TreeSummarizedExperiment(
        assays = matrix(
            ncol = nrow(object), dimnames = list(NULL, rownames(object))),
        reducedDims = list(RDA = object)
    )
    return(object)
}

# Get data for plotting
#' @importFrom scater plotReducedDim retrieveCellInfo
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
.incorporate_rda_vis <- function(
        tse, dimred, ncomponents = 2,
        colour_by = color_by, color_by = colour.by,
        colour.by = color.by, color.by = NULL,
        add.significance = TRUE, add.expl.var = TRUE, add.ellipse = TRUE,
        add.vectors = TRUE, vec.lab = NULL,
        expl.var = expl_var, expl_var = NULL,
        sep.group = "\U2014", repl.underscore = " ",
        add.centroids = FALSE, add.species = FALSE,
        # These parameters below are not used in this function. They are just
        # catched so that they are not fed into scater::plotRecucedDim.
        ellipse.alpha = 0.2, ellipse.linewidth = 0.1,
        ellipse.linetype = 1, confidence.level = 0.95,
        vec.size = 0.5, vec.color = vec.colour, vec.colour = "black",
        vec.linetype = 1, arrow.size = 0.25, min.segment.length = 5,
        label.color = label.colour, label.colour = "black", label.size = 4,
        parse.labels = TRUE, vec.text = TRUE, repel.labels = TRUE,
        position = NULL, nudge_x = NULL, nudge_y = NULL, direction = "both",
        max.overlaps = 10, check_overlap = FALSE, 
        ...){
    ###################### Input check ########################
    if( !(add.ellipse %in% c(TRUE, FALSE, "fill", "color", "colour")) ){
        stop("'add.ellipse' must be one of c(TRUE, FALSE, 'fill', ",
            "'color').", call. = FALSE)
    }
    if( !.is_a_bool(add.vectors) ){
        stop("'add.vectors must be TRUE or FALSE.", call. = FALSE)
    }
    if( !add.vectors ){
        warning("'add.vectors' is FALSE, so other arguments for vectors ",
            "and labels will be disregarded.", call. = FALSE)
    }
    if( !.is_a_bool(add.significance) ){
        stop("'add.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.expl.var) ){
        stop("'add.expl.var' must be TRUE or FALSE.", call. = FALSE)
    }
    if ( !.is_a_string(sep.group) ) {
        stop("'sep.group' must be a string specifying a separator.",
            call. = FALSE)
    }
    if ( !.is_a_string(repl.underscore) ) {
        stop("'repl.underscore' must be a string.", call. = FALSE)
    }
    if( !.is_a_bool(add.centroids) ){
        stop("'add.centroids' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(add.species) ){
        stop("'add.species' must be TRUE or FALSE.", call. = FALSE)
    }
    ###################### Input check end ####################
    # Get reducedDim
    reduced_dim <- reducedDim(tse, dimred)
    # Check that there are at least 2 coordinates
    if( ncol(reduced_dim) < 2 ){
        stop("reducedDim specified by 'dimred' must have at least 2 columns.",
            call. = FALSE)
    }
    # Only 2 dimensions are supported currently
    ncomponents <- 2
    
    # If user wants to add these, we have to get them from the actual RDA
    # object. Check if it is found
    if( add.expl.var || add.vectors || add.centroids || add.species ){
        # Check if data is available
        ind <- names(attributes(reduced_dim)) %in% c("rda", "cca", "obj")
        # If it can be found
        if( any(ind) ){
            # Add explained variance
            rda <- attributes(reduced_dim)[ind][[1]]
        } else{
            # If it cannot be found, give warning
            warning("RDA/CCA object was not found. Please compute",
                    "RDA/CCA by using addCCA or getCCA.", call. = FALSE)
            # Disable these options if the object was not found
            add.expl.var <- add.vectors  <- add.centroids <- add.species <-
                FALSE
        }
    }
    
    # If specified, get explained variance
    if( add.expl.var ){
        expl_var <- summary(rda)$concont$importance[2, ]*100
    }
    
    # Get scatter plot with plotReducedDim --> keep theme similar between
    # ordination methods
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
    # There must be at least two constrained axis to plot vectors. If there are
    # not, give warning.
    vector_data <- NULL
    if( add.vectors && ncol(rda$CCA$biplot) > 1 ){
        vector_data <- rda$CCA$biplot
        vector_data <- as.data.frame(vector_data)
        vector_data[["group"]] <- rownames(vector_data)
    } else if( add.vectors ){
        warning("Model contains only one constrained axis. Vectors cannot ",
            "be added.", call. = FALSE)
    }
    
    # Get variable names from sample metadata
    coldata <- colData(tse)
    variable_names <- colnames(coldata)
    # Check if variable names can be found metadata
    all_var_found <- FALSE
    if( !is.null(vector_data) && length(variable_names) > 0 ){
        all_var_found <- vapply(rownames(vector_data), function(x)
            vapply(variable_names, function(y) grepl(y, x), logical(1L)),
            logical(ncol(coldata)) )
        all_var_found <- all( colSums(all_var_found) == 1)
    }
    
    # Get vector labels
    if( !is.null(vector_data) ){
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
                stop("Number of labels in 'vec_lab' do not match with number ",
                    "of vectors.", call. = FALSE)
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
            warning("Significance data was not found. please compute",
                    "CCA/RDA by using addCCA or getCCA.", call. = FALSE)
        }
    }
    
    # If user wants to add site scores
    centroids <- NULL
    if( add.centroids ){
        .require_package("vegan")
        centroids <- vegan::scores(rda, display = "cn") |> as.data.frame()
        colnames(centroids) <- c("x", "y")
    }
    
    # If user wants to add species scores
    species_scores <- NULL
    if( add.species ){
        species_scores <- scores(rda, display = "species") |> as.data.frame()
        colnames(species_scores) <- c("x", "y")
    }
    # Create a list to return
    result <- list(
        plot = plot,
        ellipse_data = ellipse_data,
        vector_data = vector_data,
        centroids = centroids,
        species_scores = species_scores
    )
    return(result)
}

# Make vector labels more tidy, i.e, separate variable and group names.
# Replace also underscores with space
.tidy_vector_labels <- function(
        vector_label, coldata, sep.group, repl.underscore, ...){
    # Get variable names from sample metadata
    var_names <- colnames(coldata)
    # Loop through vector labels
    vector_label <- lapply(vector_label, FUN = function(name){
        # Get the real variable name from sample metadata
        var_name <- var_names[
            unlist(lapply(var_names, function(x) grepl(x, name))) ]
        # If the vector label includes also group name
        if( !name %in% var_names ){
            # Get the group name
            group_name <- unique( coldata[[var_name]] )[
                which(
                    paste0(var_name, unique( coldata[[var_name]] )) == name ) ]
            # Modify vector so that group is separated from variable name
            new_name <- paste0(var_name, " ", sep.group, " ", group_name)
        } else{
            new_name <- name
        }
        # Replace underscores with space
        new_name <- gsub("_", repl.underscore, new_name)
        return(new_name)
    }) |> unlist()
    return(vector_label)
}

# This function adds significance info to vector labels
.add_signif_to_vector_labels <- function(
        vector_label, var_names, signif_data, repl.underscore = " ", ...){
    # Replace underscores from significance data and variable names to match
    # labels
    rownames(signif_data) <- lapply(
        rownames(signif_data), function(x) gsub("_", repl.underscore, x)
        ) |> unlist()
    var_names <- lapply(
        var_names, function(x) gsub("_", repl.underscore, x)
        ) |> unlist()
    # Loop through vector labels
    vector_label <- lapply(vector_label, FUN = function(name){
        # Get the real variable name from sample metadata
        var_name <- var_names[
            unlist(lapply(var_names, function(x) grepl(x, name))) ]
        # Add percentage how much this variable explains, and p-value
        new_name <- expr(paste(!!name, " (",
            !!format(round(signif_data[var_name, "Explained variance"]*100, 1),
                nsmall = 1), "%, ", italic("P"), " = ",
            !!gsub("0\\.","\\.", format(round( signif_data[var_name, "Pr(>F)"],
                3), nsmall = 3)), ")"))
        return(new_name)
    }) |> unlist()
    return(vector_label)
}

# Plot based on the data
#' @importFrom ggrepel geom_text_repel geom_label_repel
.rda_plotter <- function(
        plot_data, ellipse.alpha = 0.2, ellipse.linewidth = 0.1,
        ellipse.linetype = 1, confidence.level = 0.95,
        vec.size = 0.5, vec.color = vec.colour, vec.colour = "black",
        vec.linetype = 1, arrow.size = 0.25, min.segment.length = 5,
        label.color = label.colour, label.colour = "black", label.size = 4,
        add.significance = TRUE, parse.labels = TRUE, vec.text = TRUE,
        repel.labels = TRUE, add.ellipse = TRUE,
        position = NULL, nudge_x = NULL, nudge_y = NULL, direction = "both",
        max.overlaps = 10, check_overlap = FALSE, ...){
    ###################### Input check ########################
    if( !(add.ellipse %in% c(TRUE, FALSE, "fill", "color", "colour")) ){
        stop("'add.ellipse' must be one of c(TRUE, FALSE, 'fill', ",
            "'color', 'colour').", call. = FALSE)
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
        warning("'parse.labels' was turned off because 'add.significance' ",
            "is FALSE.", call. = FALSE)
    }
    if( !.are_whole_numbers(ellipse.linetype) ){
        stop("'vec.linetype' must be a whole number.", call. = FALSE)
    }
    if ( !(is.numeric(ellipse.alpha) && ellipse.alpha > 0 &&
            ellipse.alpha < 1 ) ) {
        stop("'ellipse.alpha' must be a number between 0 and 1.",
            call. = FALSE)
    }
    if ( !(is.numeric(ellipse.linewidth) && ellipse.linewidth > 0) ) {
        stop("'ellipse.linewidth' must be a positive number.",
            call. = FALSE)
    }
    if( !(is.numeric(confidence.level) && confidence.level > 0 &&
            confidence.level < 1) ) {
        stop("'confidence.level' must be a number between 0 and 1.",
            call. = FALSE)
    }
    if ( !(is.numeric(vec.size) && vec.size > 0) ) {
        stop("'vec.size' must be a positive number.", call. = FALSE)
    }
    if ( !(is.numeric(arrow.size) && arrow.size > 0) ) {
        stop("'arrow.size' must be a positive number.", call. = FALSE)
    }
    if ( !(is.numeric(label.size) && label.size > 0) ) {
        stop("'label.size' must be a positive number.", call. = FALSE)
    }
    if ( !.is_non_empty_string(vec.color) ) {
        stop("'vec.color' must be a non-empty string specifying a colour",
            call. = FALSE)
    }
    if ( !.is_non_empty_string(label.color) ) {
        stop("'label.color' must be a non-empty string specifying a colour",
            call. = FALSE)
    }
    ###################### Input check end ####################
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
            plot <- plot + stat_ellipse(
                data = data, aes(
                    x = .data[[xvar]], y = .data[[yvar]],
                    color = .data[[colour_var]], fill = after_scale(color)),
                geom = "polygon", alpha = ellipse.alpha,
                linewidth = ellipse.linewidth, linetype = ellipse.linetype,
                level = confidence.level)
        } else if ( add.ellipse %in% c("color", "colour") ){
            plot <- plot +
                stat_ellipse(
                    data = data, aes(
                        x = .data[[xvar]], y = .data[[yvar]],
                        color = .data[[colour_var]]),
                    geom = "polygon", alpha = 0, linetype = ellipse.linetype,
                    level = confidence.level)
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
            geom_segment(
                data = data,
                aes(
                    x = 0, y = 0, xend = .data[[xvar]], yend = .data[[yvar]],
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
                box.padding = 0.25, point.padding = 1e-06, force = 1,
                force_pull = 1,
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
    
    # If user wants to add centroids
    if( !is.null(plot_data[["centroids"]]) ){
        plot <- plot +
            geom_point(
                plot_data[["centroids"]],
                mapping = aes(x = x, y = y),
                shape = 10, color = "blue")
    }
    
    # If user wants to add species scores
    if( !is.null(plot_data[["species_scores"]]) ){
        plot <- plot + geom_point(
            plot_data[["species_scores"]],
            mapping = aes(x = x, y = y),
            shape = 4, color = "red")
    }
    
    return(plot)
}
