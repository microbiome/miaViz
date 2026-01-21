
plotForest <- function(x, group.by = NULL, label.by = "CI", show.tree = TRUE){
    # Extract rowData from SE
    df <- as.data.frame(rowData(x))
    # Check label.by
    if( !is.vector(label.by) &&
        (is.null(label.by) || is.na(label.by) || label.by == "") ){
        label.by <- c()
    }
    if( !all(label.by %in% c("CI", "P-Value") | label.by %in% names(df)) ){
        stop("All terms in 'label.by' must be in rowData(x).", call. = FALSE)
    }
    # Check CI
    if( "CI" %in% label.by && !all(c("Estimate", "Lower", "Upper") %in% names(df)) ){
        stop("To show the 'CI' label, rowData(x) must include the variables ",
            "'Estimate', 'Lower' and 'Upper'.", call. = FALSE)
    }
    # Check P-Value
    if( "P-Value" %in% label.by && !"Pval" %in% names(df) ){
        stop("To show the 'P-Value' label, rowData(x) must include the ",
            "variable 'Pval'.", call. = FALSE)
    }
    # Check show.tree
    if( !is.logical(show.tree) ){
        stop("'show.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    # Plot tree
    tree.plot <- plotRowTree(
        x,
        layout = "rectangular",
        branch.length = "none",
        show.label = TRUE
    )
    # Extract tip order from tree
    tree.data <- ggplot_build(tree.plot)
    tips <- tree.data[[1]][[5]][c("y", "label")]
    tips <- tips[order(tips$y), ]
    # Order rowData based on tree tips
    tax.order <- match(tips$label, names(x))
    df <- df[tax.order, , drop = FALSE]
    # Initialise plots and widths lists
    plots <- list()
    widths <- c()
    
    if( show.tree ){
        # Store tree plot
        plots[[1]] <- plotRowTree(
            x,
            layout = "rectangular",
            branch.length = "none",
            show.label = FALSE
        )
        # Store tree width
        widths <- c(widths, 0.5)
    }
    # Convert feature names to factors
    df$Feature <- factor(rownames(df), levels = rownames(df))
    rownames(df) <- NULL
    # Find x-axis limits for forest plot
    lim <- max(abs(df$Lower), df$Upper)
    # Make forest plot
    plots[[length(plots) + 1]] <- ggplot(df, aes(x = .data$Estimate, y = .data$Feature)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "gray") +
        geom_point() +
        geom_errorbar(aes(xmin = .data$Lower, xmax = .data$Upper),
                      orientation = "y", width = 1e-2 * nrow(df)) +
        coord_cartesian(clip = "off", xlim = c(-lim, lim), ylim = c(df$Feature[1], NA)) +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              panel.grid = element_blank())
    # Store forest plot width
    widths <- c(widths, 2)
    # Construct CI label
    if( "CI" %in% label.by ){
        df$CI <- paste0(round(df$Estimate, 2), " (", round(df$Lower, 2), "—", round(df$Upper, 2), ")")
    }
    # Construct P-Value label
    if( "P-Value" %in% label.by ){
        df$`P-Value` <- paste0(round(df$Pval, 3), ifelse(df$Pval < 0.05, "*", ""))
    }
    # Initialise label plot list
    label.plots <- list()
    # Iterate over label.by terms
    for( i in seq_along(label.by) ){
        
        lab <- label.by[i]
        # Set xmax for improved positioning
        if( i == 1 ){
            xmax <- 1
        }else{
            xmax <- NA
        }
        # Make plot for current label
        label.plots[[i]] <- ggplot(df, aes(x = .data$Estimate, y = .data$Feature)) +
            geom_text(x = 0, aes(label = .data[[lab]]),
                      hjust = 0, size = 150 / nrow(df)) +
            annotate("text", x = 0, y = -1.5, label = lab, hjust = 0) +
            coord_cartesian(xlim = c(0, xmax), ylim = c(df$Feature[1], NA),
                            expand = FALSE, clip = "off") +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.background = element_blank())
    }
    # For non-empty label plot lists
    if( length(label.plots) != 0 ){
        # Wrap and store label plots
        plots[[length(plots) + 1]] <- wrap_plots(label.plots)
        # Store label plots total width
        widths <- c(widths, 1/3 * length(label.plots))
    }
    # Make final plot
    p <- wrap_plots(plots) +
        plot_layout(widths = widths)
    return(p)
}

plotForest(tse, label.by = "CI", show.tree = FALSE)
