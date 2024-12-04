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
#'   \item \code{add.vectors}: \code{Logical scalar} or \code{character vector}.
#'   If boolean, should vectors appear in the plot. If character,
#'   selects vectors that are showed. The matching is done with regular
#'   expression. (Default: \code{TRUE})
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
#'   variance appear on the coordinate axes. (Default: \code{FALSE})
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
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
setMethod("plotRDA", signature = c(x = "SingleCellExperiment"),
    function(x, dimred, ...){
        ###################### Input check ########################
        # Check dimred
        if( !dimred %in% reducedDimNames(x) ){
            stop("'dimred' must specify reducedDim.", call. = FALSE)
        }
        ###################### Input check end ####################
        # Get reducedDim
        reduced_dim <- reducedDim(x, dimred)
        # Check that there are at least 2 coordinates
        if( ncol(reduced_dim) < 2 ){
            stop("reducedDim specified by 'dimred' must have at least 2 ",
                "columns.", call. = FALSE)
        }
        # Subset by taking only constrained axes
        reduced_dim <- .subset_constrained_rda(reduced_dim)
        # Create an argument list. Only 2 dimensions are supported currently.
        args <- c(list(
            tse = x, dimred = dimred, reduced_dim = reduced_dim),
            list(...))
        args[["ncomponents"]] <- 2L
        # Get data for plotting
        plot_args <- list()
        plot_args[["ellipse_data"]] <- do.call(.get_rda_ellipse_data, args)
        plot_args[["vector_data"]] <- do.call(.get_rda_vector_data, args)
        plot_args[["centroids"]] <- do.call(.get_rda_centroids_data, args)
        plot_args[["species_scores"]] <- do.call(.get_rda_species_data, args)
        plot_args[["plot"]] <- do.call(.create_rda_baseplot, args)
        # Create a final plot
        p <- .rda_plotter(plot_args, ...)
        return(p)
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
        p <- plotRDA(x, "RDA", ...)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# Construct TreeSE from matrix to pass it to downstream functions. It is useful
# for instance if get* functios was used instead of add*.
.rda2tse <- function(object) {
    # Convert rda/cca object to TreeSE
    object <- TreeSummarizedExperiment(
        assays = SimpleList(counts = matrix(
            ncol = nrow(object), dimnames = list(NULL, rownames(object)))),
        reducedDims = list(RDA = object)
    )
    return(object)
}

# The data can include constrained and unconstrained axes. This function subsets
# the data so that it includes only constrained axes.
.subset_constrained_rda <- function(reduced_dim){
    # Get only the indices of constrained ones, i.e., first set of axes.
    # The colnames are in format, constrained_axis1, ca2, ca3..., unconstrained
    # axis1, uca2, ...
    comp_num <- as.numeric(gsub("\\D", "", colnames(reduced_dim)))
    ind <- which( cumsum(comp_num == 1) <= 1 )
    # If there were problems, it might be that the names are just arbitrary.
    # Then take all the columns.
    if( !(length(ind) > 0L && all(diff(ind) == 1L)) ){
        ind <- seq_len(ncol(reduced_dim))
    }
    # Preserve attributes
    attributes <- attributes(reduced_dim)
    attributes <- attributes[ !names(attributes) %in% c("dim", "dimnames") ]
    # Subset the data so that it includes only constrained axes
    reduced_dim <- reduced_dim[ , ind, drop = FALSE]
    if( "biplot" %in% names(attributes) ){
        attributes[["biplot"]] <- attributes[["biplot"]][ , ind, drop = FALSE]
    }
    if( "eig" %in% names(attributes) ){
        attributes[["eig"]] <- attributes[["eig"]][ind]
    }
    # Add attributes back
    attributes(reduced_dim) <- c(attributes(reduced_dim), attributes)
    return(reduced_dim)
}

# This function retrieves optional data that is used for creating an ellipses.
#' @importFrom scater retrieveCellInfo
.get_rda_ellipse_data <- function(
        tse, reduced_dim, add.ellipse = TRUE, colour_by = color_by,
        color_by = colour.by, colour.by = color.by, color.by = NULL, ...){
    #
    if( !(add.ellipse %in% c(TRUE, FALSE, "fill", "color", "colour") &&
            length(add.ellipse) == 1L ) ){
        stop("'add.ellipse' must be one of c(TRUE, FALSE, 'fill', ",
            "'color').", call. = FALSE)
    }
    if( !(is.null(colour_by) || .is_a_string(colour_by) &&
            colour_by %in% colnames(colData(tse)) ) ){
        stop("'colour_by' must be NULL or name of column from colData(x).",
            call. = FALSE)
    }
    #
    ellipse_data <- NULL
    if( add.ellipse != FALSE && !is.null(colour_by) ){
        # Ellipse data is the same ordination data
        ellipse_data <- as.data.frame(reduced_dim)
        # Add sample metadata from colData
        ellipse_data[[colour_by]] <- retrieveCellInfo(tse, colour_by)[["value"]]
        attributes(ellipse_data)[["colour_by"]] <- colour_by
    }
    return(ellipse_data)
}

# This function retrieves data for creating vectors. Moreover, it wrangles the
# vector data and controls what information is added to vector text or labels.
.get_rda_vector_data <- function(
        tse, reduced_dim, add.vectors = TRUE, add.significance = TRUE,
        vec.lab = NULL, sep.group = "\U2014", repl.underscore = " ",
        ignore.case = FALSE, ...){
    #
    if( !( .is_a_bool(add.vectors) || is.character(add.vectors)) ){
        stop("'add.vectors must be TRUE or FALSE or character vector.",
            call. = FALSE)
    }
    if( !.is_a_bool(add.significance) ){
        stop("'add.significance' must be TRUE or FALSE.", call. = FALSE)
    }
    if ( !.is_a_string(sep.group) ) {
        stop("'sep.group' must be a string specifying a separator.",
            call. = FALSE)
    }
    if ( !.is_a_string(repl.underscore) ) {
        stop("'repl.underscore' must be a string.", call. = FALSE)
    }
    if( !.is_a_bool(ignore.case) ){
        stop("'ignore.case' must be TRUE or FALSE.", call. = FALSE)
    }
    if( add.significance && (.is_a_bool(add.vectors) && !add.vectors) ){
        # If it cannot be found, give warning
        warning("'add.vectors' is FALSE, so other arguments for vectors and ",
                "labels will be disregarded.", call. = FALSE)
    }
    #
    # There must be at least two constrained axis to plot vectors. If there are
    # not, give warning.
    if( ncol(reduced_dim) <= 1 ){
        add.vectors <- FALSE
        warning("Model contains only one constrained axis. Vectors cannot ",
                "be added.", call. = FALSE)
    }
    # Get vector data, i.e, biplot
    vector_data <- if(!(.is_a_bool(add.vectors) && !add.vectors))
        .get_rda_attribute(reduced_dim, "biplot")
    if( !is.null(vector_data) ){
        vector_data <- as.data.frame(vector_data)
        vector_data[["group"]] <- rownames(vector_data)
        vector_data[["vector_label"]] <- rownames(vector_data)
        # Subset vectors; show only specified ones, if add.vectors specifies the
        # covariate names. The matching is done with regular expression so that
        # user can specify for instance the whole covariate easily (covariate
        # names and group values are merged in biplot).
        if( is.character(add.vectors) ){
            add.vectors <- paste0(add.vectors, collapse = "|")
            keep <- vapply(rownames(vector_data), function(x)
                grepl(add.vectors, x, perl = TRUE, ignore.case = ignore.case),
                logical(1L))
            vector_data <- vector_data[keep, ]
        }
        # If all vectors were removed, give NULL
        if( nrow(vector_data) == 0L ){
            vector_data <- NULL
        }
    }
    
    # Get sample metadata. Check if all biplot covariate names can be found
    # from sample metadata. As biplot have merged names, we have to use sample
    # metadata later to make the vector labels tidier.
    coldata <- colData(tse)
    variable_names <- colnames(coldata)
    all_var_found <- vapply(rownames(vector_data), function(x)
        vapply(variable_names, function(y) grepl(y, x), logical(1L)),
        logical(ncol(coldata)) )
    all_var_found <- all( colSums(all_var_found) == 1)
    
    # Make the vector labels tidier. For instance, covriate name and value
    # are separated. This applies only when labels were not provided by user.
    if( !is.null(vector_data) && is.null(vec.lab) && all_var_found ){
        # Make labels more tidy
        vector_data[["vector_label"]] <- .tidy_vector_labels(
            vector_data[["vector_label"]], coldata, sep.group = sep.group,
            repl.underscore = repl.underscore)
    }
    # Add vector labels provodied by user
    if( !is.null(vector_data) && !is.null(vec.lab) ){
        # Check that user-provided labels are correct length
        if( length(vec.lab) != nrow(vector_data) ){
            stop("Number of labels in 'vec_lab' do not match with number ",
                "of vectors.", call. = FALSE)
        }
        # If they are, add labels to data
        vector_data[["vector_label"]] <- vec.lab
    }
    
    # Add significance information to the labels
    signif_data <- if( add.significance && !is.null(vector_data) &&
        all_var_found ) .get_rda_attribute(reduced_dim, "significance")
    if( !is.null(signif_data) ){
        signif_data <- signif_data[[1L]] |> as.data.frame()
        # Add significance to vector labels
        vector_data[["vector_label"]] <- .add_signif_to_vector_labels(
            vector_data[["vector_label"]], variable_names, signif_data,
            repl.underscore)
    }
    if( add.significance && is.null(signif_data) &&
            !(.is_a_bool(add.vectors) && !add.vectors) ){
        # If it cannot be found, give warning
        warning("Significance data was not found. please compute",
                "CCA/RDA by using add* function.", call. = FALSE)
    }
    
    return(vector_data)
}

# Make vector labels more tidy, i.e, separate variable and group names. We need
# colData for this because in the biplot data, covariate name and values are
# just combined, i.e., we do not know if "group" is the group name or is it
# "groupName" when the name in biplot is "groupNameValue". 
# Replace also underscores with space.
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
        new_name <- expr(
            paste(!!name, " (",
                !!format(
                    round(signif_data[var_name, "Explained variance"]*100, 1),
                    nsmall = 1), "%, ", italic("P"), " = ",
                !!gsub("0\\.","\\.", format(
                    round(signif_data[var_name, "Pr(>F)"], 3),
                    nsmall = 3)), ")"))
        return(new_name)
    }) |> unlist()
    return(vector_label)
}

# This functions returns optional centroids for plotting.
.get_rda_centroids_data <- function(
        reduced_dim, add.centroids = FALSE, ncomponents = 2L, ...){
    #
    if( !.is_a_bool(add.centroids) ){
        stop("'add.centroids' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_an_integer(ncomponents) ){
        stop("'ncomponents' must be an integer.", call. = FALSE)
    }
    #
    res <- if(add.centroids) .get_rda_attribute(reduced_dim, "centroids")
    if( !is.null(res) ){
        res <- res[, seq_len(ncomponents), drop = FALSE] |> as.data.frame()
        colnames(res) <- c("x", "y")
    }
    return(res)
}

# This functions returns optional species scores for plotting.
.get_rda_species_data <- function(
        reduced_dim, add.species = FALSE, ncomponents = 2L, ...){
    if( !.is_a_bool(add.species) ){
        stop("'add.species' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_an_integer(ncomponents) ){
        stop("'ncomponents' must be an integer.", call. = FALSE)
    }
    #
    res <- if(add.species) .get_rda_attribute(reduced_dim, "species")
    if( add.species ){
        res <- res[, seq_len(ncomponents), drop = FALSE] |> as.data.frame()
        colnames(res) <- c("x", "y")
    }
    return(res)
}

# This function is used to fetch specified datatype from attributes if it
# exists.
.get_rda_attribute <- function(reduced_dim, attr_names){
    res <- NULL
    attr_values <- attributes(reduced_dim)
    ind <- names(attr_values) %in% attr_names
    if( any(ind) ){
        res <- attr_values[ind][[1L]]
    }
    return(res)
}

# This function utilizes scater::plotReducedDim to create "baseplot". To where
# we can build the the plot. The idea is that the theme is similar in all
# ordination plots.
#' @importFrom scater plotReducedDim
.create_rda_baseplot <- function(
        tse, dimred, reduced_dim, ncomponents = 2L,
        add.expl.var = FALSE, expl.var = expl_var, expl_var = NULL,
        colour_by = color_by, color_by = colour.by,
        colour.by = color.by, color.by = NULL, ...){
    #
    if( !.is_a_bool(add.expl.var) ){
        stop("'add.expl.var' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_an_integer(ncomponents) ){
        stop("'ncomponents' must be an integer.", call. = FALSE)
    }
    if( !( is.null(expl.var) || (is.numeric(expl.var) &&
            length(expl.var) == ncomponents )) ){
        stop("'expl.var' must be numeric vector with length ", ncomponents,
            ".", call. = FALSE)
    }
    if( !(is.null(colour_by) || .is_a_string(colour_by) &&
            colour_by %in% colnames(colData(tse)) ) ){
        stop("'colour_by' must be NULL or name of column from colData(x).",
            call. = FALSE)
    }
    #
    # If specified, get explained variance
    if( add.expl.var && !is.null(expl.var) ){
        eigen_vals <- attr(reduced_dim, "eig")
        # Convert to explained variance and take only first two components
        expl_var <- eigen_vals / sum(eigen_vals)
        expl_var <- expl_var[seq_len(ncomponents)]*100
        expl_var <- summary(rda)$concont$importance*100
    }
    # Create argument list
    args <- c(list(object = tse, dimred = dimred, ncomponents = ncomponents,
        colour_by = colour_by, percentVar = expl_var), list(...))
    # Remove additional arguments since plotReducedDim fails if we feed
    # values that are not recognized
    remove <- names(args) %in% c(
        "add.significance", "add.expl.var", "add.ellipse", "add.vectors",
        "vec.lab", "sep.group", "repl.underscore", "add.centroids",
        "add.species", "ellipse.alpha", "ellipse.linewidth", "ellipse.linetype",
        "confidence.level", "vec.size", "vec.color", "vec.colour",
        "vec.linetype", "arrow.size", "min.segment.length", "label.color",
        "label.colour", "label.size", "parse.labels", "vec.text",
        "repel.labels", "position", "nudge_x", "nudge_y", "direction",
        "max.overlaps", "check_overlap", "ignore.case")
    args <- args[ !remove ]
    # Get scatter plot with plotReducedDim --> keep theme similar between
    # ordination methods
    p <- do.call(plotReducedDim, args)
    return(p)
}

# This function is used to create the plot.
.rda_plotter <- function(plot_data, ...){
    # Get the scatter plot
    plot <- plot_data[["plot"]]
    # Add ellipse
    plot <- .rda_plotter_ellipse(plot, plot_data, ...)
    # Add vectors
    plot <- .rda_plotter_vector(plot, plot_data, ...)
    # Add centroids
    plot <- .rda_plotter_centroids_or_species(plot, plot_data, "centroids")
    # Add species
    plot <- .rda_plotter_centroids_or_species(plot, plot_data, "species_scores")
    return(plot)
}

# This function adds ellipse visualization.
.rda_plotter_ellipse <- function(
        plot, plot_data, add.ellipse = TRUE, ellipse.alpha = 0.2,
        ellipse.linewidth = 0.1, ellipse.linetype = 1, confidence.level = 0.95,
        ...){
    #
    if( !(add.ellipse %in% c(TRUE, FALSE, "fill", "color", "colour") &&
            length(add.ellipse) == 1L ) ){
        stop("'add.ellipse' must be one of c(TRUE, FALSE, 'fill', ",
            "'color').", call. = FALSE)
    }
    if( !.are_whole_numbers(ellipse.linetype) ){
        stop("'vec.linetype' must be a whole number.", call. = FALSE)
    }
    if ( !(is.numeric(ellipse.alpha) && ellipse.alpha > 0 &&
            ellipse.alpha < 1 ) ) {
        stop("'ellipse.alpha' must be a number between 0 and 1.", call. = FALSE)
    }
    if ( !(is.numeric(ellipse.linewidth) && ellipse.linewidth > 0) ) {
        stop("'ellipse.linewidth' must be a positive number.", call. = FALSE)
    }
    if( !(is.numeric(confidence.level) && confidence.level > 0 &&
            confidence.level < 1) ) {
        stop("'confidence.level' must be a number between 0 and 1.",
            call. = FALSE)
    }
    #
    data <- plot_data[["ellipse_data"]]
    if( !is.null(data) ){
        xvar <- colnames(data)[[1]]
        yvar <- colnames(data)[[2]]
        colour_var <- attributes(data)[["colour_by"]]
        # Add ellipses to plot (fill or colour the edge)
        fill <- add.ellipse %in% c(TRUE, "fill")
        plot <- plot + stat_ellipse(
            data = data,
            mapping = aes(
                x = .data[[xvar]], y = .data[[yvar]],
                color = .data[[colour_var]], fill = after_scale(color)),
            geom = "polygon",
            linewidth = ellipse.linewidth,
            linetype = ellipse.linetype,
            level = confidence.level,
            alpha = if(fill) ellipse.alpha else 0
            )
    }
    return(plot)
}

# This function adds vector and text layer to the plot.
#' @importFrom ggrepel geom_text_repel geom_label_repel
.rda_plotter_vector <- function(
        plot, plot_data, vec.size = 0.5, arrow.size = 0.25, label.size = 4,
        vec.color = vec.colour, vec.colour = "black",
        label.color = label.colour, label.colour = "black",
        vec.text = TRUE, repel.labels = TRUE, parse.labels = TRUE,
        add.significance = TRUE, vec.linetype = 1, min.segment.length = 5,
        position = NULL, nudge_x = NULL, nudge_y = NULL, direction = "both",
        max.overlaps = 10, check_overlap = FALSE, ...
        ){
    #
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
    #
    data <- plot_data[["vector_data"]]
    if( !is.null(data) ){
        xvar <- colnames(data)[[1]]
        yvar <- colnames(data)[[2]]
        # Add vectors
        plot <- plot +
            geom_segment(data = data, aes(
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
        # Based on desired outcome, get correct argument set
        if( repel.labels ){
            label_args <- c(
                label_args, min.segment.length = min.segment.length,
                box.padding = 0.25, point.padding = 1e-06, force = 1,
                force_pull = 1,
                max.time = 0.5, max.iter = 10000, max.overlaps = max.overlaps,
                direction = direction, seed = NA, verbose = FALSE
            )
        } else if( !repel.labels && vec.text ){
            label_args <- c(label_args, check_overlap = check_overlap)
        }
        # Choose right function and call it
        FUN <- if( repel.labels && vec.text ) geom_text_repel
            else if( repel.labels && !vec.text ) geom_label_repel
            else if( !repel.labels && vec.text ) geom_text
            else geom_label
        plot <- plot + do.call(FUN, label_args)
    }
    return(plot)
}

# This function adds centroids or species layer to the plot.
.rda_plotter_centroids_or_species <- function(plot, plot_data, type){
    data <- plot_data[[type]]
    if( !is.null(data) ){
        plot <- plot + geom_point(
            data,
            mapping = aes(x = x, y = y),
            shape = if(type == "centroids") 10L else 4L,
            color = if(type == "centroids") "blue" else "red",
            )
    }
    return(plot)
}
