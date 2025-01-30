#' @name
#' plotMediation
#'
#' @title
#' Visualize mediation
#'
#' @description
#' \code{plotMediation()} generates a heatmap from the results of mediation
#' analysis produced with \code{mia:getMediation()} or
#' \code{mia:addMediation()}. It displays effect size and significance of the
#' Actual Causal Mediation Effect (ACME) and the Actual Direct Effect (ADE) for
#' each mediator included in the object \code{x}.
#'
#' @details
#' \code{plotMediation} creates a heatmap starting from the
#' output of the \code{\link[mia:getMediation]{mediation}} functions, which are
#' mia wrappers for the basic \code{\link[mediation:mediate]{mediate}} function.
#' Either a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' or a data.frame object is supported as input. When the input is a
#' SummarizedExperiment, this should contain the output of addMediation
#' in the metadata slot and the argument \code{name} needs to be defined.
#' When the input is a data.frame, this should be returned as output from
#' getMediation.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object or a \code{data.frame}, returned as output from
#' \code{\link[mia:getMediation]{addMediation}} or
#' \code{\link[mia:getMediation]{getMediation}}, respectively.
#'
#' @param name \code{Character scalar} value defining which mediation data
#' to use. (Default: \code{"mediation"})
#'
#' @param layout \code{Character scalar} Determines the layout of plot. Must be
#' either "heatmap" or "forest". (Default: \code{"heatmap"})
#'
#' @param ... Additional parameters for plotting.
#' \itemize{
#'   \item \code{add.significance}: \code{Character scalar}. Controls how
#'   p-values are displayed in the heatmap. Options include \code{"symbol"}
#'   (p-values displayed as symbols), \code{"numeric"} (p-values displayed as
#'   numeric values) and \code{NULL} (p-values not displayed).
#'   (Default: \code{"symbol"})
#' }
#'
#' @examples
#' \dontrun{
#' library(mia)
#' library(scater)
#'
#' # Load dataset
#' data(hitchip1006, package = "miaTime")
#' tse <- hitchip1006
#'
#' # Agglomerate features by family (merely to speed up execution)
#' tse <- agglomerateByRank(tse, rank = "Phylum")
#' # Convert BMI variable to numeric
#' tse[["bmi_group"]] <- as.numeric(tse[["bmi_group"]])
#'
#' # Apply clr transformation to counts assay
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#'
#' # Analyse mediated effect of nationality on BMI via clr-transformed features
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' tse <- addMediation(
#'     tse, name = "assay_mediation",
#'     outcome = "bmi_group",
#'     treatment = "nationality",
#'     assay.type = "clr",
#'     covariates = c("sex", "age"),
#'     treat.value = "Scandinavia",
#'     control.value = "CentralEurope",
#'     boot = TRUE, sims = 100,
#'     p.adj.method = "fdr"
#'     )
#'
#' # Visualise results as heatmap
#' plotMediation(tse, "assay_mediation")
#'
#' # Visualise results as forest plot
#' plotMediation(tse, "assay_mediation", layout = "forest")
#'}
#'
NULL

#' @rdname plotMediation
#' @export
#' @importFrom S4Vectors metadata
setMethod("plotMediation", signature = c(x = "SummarizedExperiment"),
    function(x, name = "mediation", ...){
        # Check metadata name
        if( !(.is_a_string(name) && name %in% names(metadata(x))) ){
            stop("'name' must be in names(metadata(x)).", call. = FALSE)
        }
        #
        df <- metadata(x)[[name]]
        p <- plotMediation(df, ...)
        return(p)
    }
)

#' @rdname plotMediation
#' @export
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr rename_with
setMethod("plotMediation", signature = c(x = "data.frame"),
    function(x, layout = "heatmap", ...){
        # Check layout
        if( !(.is_a_string(layout) && layout %in% c("heatmap", "forest")) ){
            stop("'layout' must be either 'heatmap' or 'forest'.",
                call. = FALSE)
        }
        # Wrangle meditation data into correct format
        df <- .get_mediation_data(x, ...)
        # Create the plot
        FUN <- switch(layout,
            "heatmap" = .plot_mediation_heatmap,
            "forest" = .plot_mediation_forest
            )
        p <- FUN(df)
        return(p)
    }
)

################################ HELP FUNCTIONS ################################

# This function converts mediation table into long format, ready for plotting.
#' @importFrom dplyr rename_with mutate case_when
#' @importFrom tidyr pivot_longer pivot_wider
.get_mediation_data <- function(df, add.significance = "symbol", ...){
    # To disable "no visible binding for global variable" message in cmdcheck
    acme <- ade <- total <- mediator <- Metric <- value <- NULL

    # Check add.significance
    if( !(is.null(add.significance) || (.is_a_string(add.significance))
            && add.significance %in% c("symbol", "numeric") )){
        stop("'add.significance' must be NULL, 'symbol' or 'numeric'.",
            call. = FALSE)
    }
    # Check that the right columns can be found, i.e., mediation is calculated
    # with mia.
    cols <- c("acme", "ade", "total")
    cols <- paste0(
        cols,
        rep(c("", "_lower", "_upper", "_pval", "_padj"), each = length(cols)))
    if( !all(cols %in% colnames(df)) ){
        stop("Mediation results seems to be in wrong format. Please calculate ",
            "the results with mia::getMediation() or mia::addMediation().",
            call. = FALSE)
    }
    # Put the data in correct format. Estimate, condifednce intervals and
    # p-values goes to separate columns. ACME, ADE and total esimates are in
    # long format.
    df <- df |>
        # Rename to add_estimate suffix which is next used to categorize values
        rename_with(~ paste0(.x, "_estimate"), c(acme, ade, total)) |>
        # Convert into long format. Get metric type from suffix
        pivot_longer(
            cols = -mediator, names_to = c("Condition", "Metric"),
            names_sep = "_", names_prefix = "hi") |>
        # Estimates, confidence intervals and p-values into own columns
        pivot_wider(names_from = Metric, values_from = value)
    # Convert to factor to preserve the order when plotting
    df[["Condition"]] <- toupper( df[["Condition"]] )
    df[["Condition"]] <- factor(
        df[["Condition"]], levels = rev(unique(df[["Condition"]])))
    # Sort based on padj
    df[["mediator"]] <- factor(
        df[["mediator"]],
        levels = unique(df[["mediator"]][order(df[["padj"]])]))
    # Control p-value layout
    if( is.null(add.significance) ){
        # Do not add significance, i.e., add empty label
        df[["padj"]] <- rep("", nrow(df))
    } else if( add.significance == "symbol" ){
        # Add p-values as symbol
        df <- df |>
            mutate(padj = case_when(
                padj < 0.001 ~ "***",
                padj < 0.010 ~ "**",
                padj < 0.050 ~ "*",
                TRUE ~ ""
            ))
    } else{
        # If p-values are displayed as raw values, round them and convert to
        # character
        df[["padj"]] <- as.character( round(df[["padj"]], 3L) )
    }
    return(df)
}

# This function creates a heatmap from mediation results
.plot_mediation_heatmap <- function(df){
    .require_package("shadowtext")
    # To disable "no visible binding for global variable" message in cmdcheck
    Condition <- mediator <- estimate <- padj <- NULL

    # For heatmap color scale
    effect_max <- max(abs(df[["estimate"]]))
    #
    p <- ggplot(df, aes(x = Condition, y = mediator, fill = estimate)) +
        # Add heatmap tiles
        geom_tile(
            mapping = aes(fill = estimate),
            color = "white", lwd = 1.5, linetype = 1) +
        # Adjust color scale
        scale_fill_gradientn(
            colours = c("blue", "white", "red"),
            limits = c(-effect_max, effect_max), name = "Effect") +
        # Add p-values
        shadowtext::geom_shadowtext(aes(label = padj)) +
        # Modify x and y axis titles
        labs(x = "", y = "")
    return(p)
}

# This function creates a forest plot from mediation results.
.plot_mediation_forest <- function(df){
    # To disable "no visible binding for global variable" message in cmdcheck
    estimate <- Condition <- lower <- upper <- NULL

    #
    p <- ggplot(df, aes(x = estimate, y = Condition)) +
        # Add points for estimates
        geom_point(size = 3) +
        # Add error bars for confidence intervals
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
        # Add a line to 0 to denote zero effect
        geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
        facet_wrap(. ~ mediator) +
        # Modify x and y axis titles
        labs(x = "Effect", y = "")
    return(p)
}
