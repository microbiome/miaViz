library(ggord)
library(ggplot2)

# import(ggord)
# import(ggplot2)

plotRDA <- function(tse, colour_by = NULL, color_by = colour_by) {
  
  # Extract RDA object from reducedDims of tse
  rda <- reducedDim(tse, "RDA")
  
  # Store permanova test results into variable
  signif_info <- attr(rda, "significance")$permanova
  
  # Get dimensions
  rda <- attr(rda, "rda")
  
  # Adjust names
  # Get labels of vectors
  vec_lab_old <- rownames(rda$CCA$biplot)
  
  # Loop through vector labels
  vec_lab <- sapply(vec_lab_old, FUN = .name_var, tse = tse, signif_info = signif_info)
  
  # Add names
  names(vec_lab) <- vec_lab_old
  
  # Create labels for axis
  xlab <- paste0("RDA1 (", format(round( rda$CCA$eig[[1]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
  ylab <- paste0("RDA2 (", format(round( rda$CCA$eig[[2]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
  
  # If colour var is stated, ggord grouping takes the correspondig column from
  # coldata, else ggord grouping is set to NULL
  if ( !is.null(colour_by) ) {
    grp_in <- colData(tse)[[colour_by]]
  } else {
    grp_in <- NULL
  }
  
  # Create a plot        
  plot <- ggord(rda,
                grp_in = grp_in,
                vec_lab = vec_lab,
                alpha = 0.5,
                size = 4, addsize = -4,
                #ext= 0.7, 
                #coord_fix = FALSE,
                txt = 3.5, repel = TRUE) + 
    # Adjust titles
    theme( axis.title = element_text(size = 10) )
  
  if ( !is.null(colour_by) ) {
    
    # Adjust labels
    plot <- plot +
      guides(colour = guide_legend(colour_by),
      fill = guide_legend(colour_by),
      group = guide_legend(colour_by),
      shape = guide_legend(colour_by),
      x = guide_axis(xlab),
      y = guide_axis(ylab))
    
  }
  
  return(plot)

}


.name_var <- function(name, tse, signif_info){
  
  # Retrieve model variable names
  variable_names <- rownames(signif_info)[2:(nrow(signif_info) - 1)]
  # Get the variable name
  variable_name <- grep(paste(variable_names, collapse = "|"), name, value = TRUE)
  # If the vector label includes also group name
  if( !any(name %in% variable_names) ){
    
    # Get the group names
    group_name <- unique( colData(tse)[[variable_name]] )[ which( paste0(variable_name, unique( colData(tse)[[variable_name]] )) == name ) ]
    
    # Modify vector so that group is separated from variable name
    new_name <- paste0(variable_name, " \U2012 ", group_name)
    
  } else {
    new_name <- name
  }
  
  # Add percentage how much this variable explains, and p-value
  new_name <- expr(paste(!!new_name, " (", 
                         !!format(round( signif_info[variable_name, "Explained variance"]*100, 1), nsmall = 1), 
                         "%, ",italic("P"), " = ", 
                         !!gsub("0\\.","\\.", format(round( signif_info[variable_name, "Pr(>F)"], 3), 
                                                     nsmall = 3)), ")"))
  
  return(new_name)

}

# rda_plot <- plotRDA(rda_tse, colour_by = "ClinicalStatus")
