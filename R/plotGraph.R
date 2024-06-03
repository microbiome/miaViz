#' Plotting igraph objects with information from a \code{SummarizedExperiment}
#' 
#' \code{plotGraph} plots an \code{igraph} object with additional information
#' matched from a \code{SummarizedExperiment} object for the nodes only.
#' Information on the edges have to provided manually.
#'
#' @param x,y a graph object and a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   object or just a 
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'   For the latter object a graph object must be stored in \code{metadata(x)$name}.
#'   
#' @param name If \code{x} is a 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' the key for subsetting the \code{metadata(x)} to a graph object.
#' 
#' @param show.label \code{logical} (scalar), \code{integer} or \code{character}
#'   vector. If a \code{logical} scalar is given, should tip labels be plotted
#'   or if a logical vector is provided, which labels should be shown? If an
#'   \code{integer} or \code{character} vector is provided, it will be converted
#'   to a logical vector. The \code{integer} values must be in the range of 1
#'   and number of nodes, whereas the values of a \code{character} vector must
#'   match values of a \code{label} or \code{name} column in the node data. In
#'   case of a \code{character} vector only values corresponding to actual
#'   labels will be plotted and if no labels are provided no labels will be
#'   shown. (default: \code{show.label = FALSE})
#'   
#' @param show_label Deprecated. Use \code{show.label} instead.
#' 
#' @param add_legend logical scalar. Should legends be plotted? 
#'   (default: \code{add_legend = TRUE})
#'   
#' @param layout layout for the plotted graph. See 
#'   \code{\link[ggraph:ggraph]{ggraph}} for details. (default: 
#'   \code{layout = "kk"})
#'   
#' @param edge.type type of edge plotted on the graph. See 
#'   \code{\link[ggraph:geom_edge_fan]{geom_edge_fan}} for details and other 
#'   available geoms. (default: 
#'   \code{edge.type = "fan"})
#' 
#' @param edge_type Deprecated. Use \code{edge.type} instead.
#'   
#' @param edge.colour.by Specification of a edge metadata field to use for 
#'   setting colours of the edges.
#'   
#' @param edge_colour_by Deprecated. Use \code{edge.colour.by} instead.
#'   
#' @param edge.width.by Specification of a edge metadata field to use for 
#'   setting width of the edges.
#'   
#' @param edge_width_by Deprecated. Use \code{edge.width.by} instead.
#' 
#' @param colour_by Specification of a column metadata field or a feature to
#'   colour graph nodes by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param shape.by Specification of a column metadata field or a feature to
#'   shape graph nodes by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param shape_by Deprecated. Use \code{shape.by} instead.
#'   
#' @param size.by Specification of a column metadata field or a feature to
#'   size graph nodes by, see the by argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible 
#'   values.
#'   
#' @param size_by Deprecated. Use \code{size.by} instead.
#'   
#' @param assay.type A string or integer scalar specifying which assay to
#'   obtain expression values from, for use in point aesthetics - see the 
#'   \code{exprs_values} argument in 
#'   \code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}}.
#'   
#' @param other.fields Additional fields to include in the node information
#'   without plotting them.
#'   
#' @param other_fields Deprecated. Use \code{other.fields} instead. 
#'   
#' @param ... additional arguments for plotting. See 
#'   \code{\link{mia-plot-args}} for more details i.e. call \code{help("mia-plot-args")}
#' 
#' @details:
#' Internally \code{tidygraph} and \code{ggraph} are used. Therefore, all 
#' graph types which can be converted by \code{tidygraph::as_tbl_graph} can
#' be used.
#' 
#' @return a \code{\link{ggtree}} plot
#' 
#' @name plotGraph
#' 
#' @examples
#' \donttest{
#' # data setup
#' library(mia)
#' data(GlobalPatterns)
#' data(col_graph)
#' data(row_graph)
#' data(row_graph_order)
#' metadata(GlobalPatterns)$col_graph <- col_graph
#' 
#' genus <- agglomerateByRank(GlobalPatterns,"Genus",na.rm=TRUE)
#' metadata(genus)$row_graph <- row_graph
#' order <- agglomerateByRank(genus,"Order",na.rm=TRUE)
#' metadata(order)$row_graph <- row_graph_order
#' 
#' # plot a graph independently
#' plotColGraph(col_graph,
#'              genus,
#'              colour_by = "SampleType",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight",
#'              show.label = TRUE)
#' 
#' # plot the graph stored in the object
#' plotColGraph(genus,
#'              name = "col_graph",
#'              colour_by = "SampleType",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight")
#'              
#' 
#' # plot a graph independently
#' plotRowGraph(row_graph,
#'              genus,
#'              colour_by = "Kingdom",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight")
#' 
#' # plot the graph stored in the object
#' plotRowGraph(genus,
#'              name = "row_graph",
#'              colour_by = "Phylum",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight")
#'
#'                            
#' # plot a graph independently
#' plotRowGraph(row_graph_order,
#'              order,
#'              colour_by = "Kingdom",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight")
#' 
#' # plot the graph stored in the object and include some labels
#' plotRowGraph(order,
#'              name = "row_graph",
#'              colour_by = "Phylum",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight", 
#'              show.label = c("Sulfolobales","Spirochaetales",
#'                             "Verrucomicrobiales"))
#'                             
#' # labels can also be included via selecting specific rownames of x/y
#' plotRowGraph(order,
#'              name = "row_graph",
#'              colour_by = "Phylum",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight", 
#'              show.label = c(1,10,50))
#'              
#' # labels can also be included via a logical vector, which has the same length
#' # as nodes are present
#' label_select <- rep(FALSE,nrow(order))
#' label_select[c(1,10,50)] <-  TRUE
#' plotRowGraph(order,
#'              name = "row_graph",
#'              colour_by = "Phylum",
#'              edge.colour.by = "weight",
#'              edge.width.by = "weight",
#'              show.label = label_select)
#' }
NULL

#' @rdname plotGraph
#' @export
setGeneric("plotColGraph", signature = c("x","y"),
           function(x, y, ...) standardGeneric("plotColGraph"))

#' @rdname plotGraph
#' @export
setGeneric("plotRowGraph", signature = c("x","y"),
           function(x, y, ...) standardGeneric("plotRowGraph"))

.check_graph_plot_switches <- function(show.label, add_legend){
    if(!.is_a_bool(show.label)){
        if( (!is.logical(show.label) && !is.character(show.label) && 
            !is.numeric(show.label)) ||
            is.null(show.label)){
            stop("'show.label' must be either TRUE or FALSE or logical, ",
                 "integer or character ",
                 "vector. Character alues should match the label of the graph.",
                 call. = FALSE)
        }
    }
    if(!.is_a_bool(add_legend)){
        stop("'add_legend' must be either TRUE or FALSE.", call. = FALSE)
    }
}

.norm_layout_edge_type <- function(layout, edge.type){
    edge.type <- match.arg(edge.type[1L], c("fan","link","arc","parallel"))
    return(list(layout = layout,
                edge.type = edge.type))
}

#' @rdname plotGraph
#' @importFrom tidygraph as_tbl_graph
#' @export
setMethod("plotColGraph",
    signature = c(x = "ANY",y = "SummarizedExperiment"),
    function(x, y,
             show.label = show_label,
             show_label = FALSE,
             add_legend = TRUE,
             layout = "kk",
             edge.type = c("fan","link","arc","parallel"),
             edge.colour.by = edge_colour_by,
             edge_colour_by = NULL,
             edge.width.by = edge_width_by,
             edge_width_by = NULL,
             colour_by = NULL,
             shape.by = shape_by,
             shape_by = NULL,
             size.by = size_by,
             size_by = NULL,
             assay.type = "counts",
             other.fields = other_fields,
             other_fields = list(),
             ...){
        .plot_row_column_graph(x = x, y = y,
                              show.label = show.label,
                              add_legend = add_legend,
                              layout = layout,
                              edge.type = edge.type,
                              edge.colour.by = edge.colour.by,
                              edge.width.by = edge.width.by,
                              colour_by = colour_by,
                              shape.by = shape.by,
                              size.by = size.by,
                              assay.type = assay.type,
                              other.fields = other.fields,
                              type = "column",
                              ...)
    }
)

#' @rdname plotGraph
#' @importFrom S4Vectors metadata
#' @export
setMethod("plotColGraph",
    signature = c(x = "SummarizedExperiment", y = "missing"),
    function(x, y, name = "graph", ...){
        graph <- metadata(x)[[name]]
        if(is.null(graph)){
            stop("No data found in metadata for key '",name,"'", call. = FALSE)
        }
        plotColGraph(graph, x, ...)
    }
)

#' @rdname plotGraph
#' @importFrom tidygraph as_tbl_graph
#' @export
setMethod("plotRowGraph",
    signature = c(x = "ANY",y = "SummarizedExperiment"),
    function(x, y,
             show.label = show_label,
             show_label = FALSE,
             add_legend = TRUE,
             layout = "kk",
             edge.type = c("fan","link","arc","parallel"),
             edge.colour.by = edge_colour_by,
             edge_colour_by = NULL,
             edge.width.by = edge_width_by,
             edge_width_by = NULL,
             colour_by = NULL,
             shape.by = shape_by,
             shape_by = NULL,
             size.by = NULL,
             assay.type = "counts",
             other.fields = other_fields,
             other_fields = list(),
             ...){
        .plot_row_column_graph(x = x, y = y,
                               show.label = show.label,
                               add_legend = add_legend,
                               layout = layout,
                               edge.type = edge.type,
                               edge.colour.by = edge.colour.by,
                               edge.width.by = edge.width.by,
                               colour_by = colour_by,
                               shape.by = shape.by,
                               size.by = size.by,
                               assay.type = assay.type,
                               other.fields = other.fields,
                               type = "row",
                               ...)
    }
)

#' @rdname plotGraph
#' @importFrom S4Vectors metadata
#' @export
setMethod("plotRowGraph",
    signature = c(x = "SummarizedExperiment", y = "missing"),
    function(x, y, name = "graph",...){
        graph <- metadata(x)[[name]]
        if(is.null(graph)){
            stop("No data found in metadata for key '",name,"'", call. = FALSE)
        }
        plotRowGraph(graph, x, ...)
    }
)

.plot_row_column_graph <- function(x, y,
                                   show.label = FALSE,
                                   add_legend = TRUE,
                                   layout = "kk",
                                   edge.type = c("fan","link","arc","parallel"),
                                   edge.colour.by = NULL,
                                   edge.width.by = NULL,
                                   colour_by = NULL,
                                   shape.by = NULL,
                                   size.by = NULL,
                                   assay.type = "counts",
                                   other.fields = other_fields,
                                   other_fields = list(),
                                   type = c("row","column"),
                                   ...){
    type <- match.arg(type)
    # input check
    .check_graph_plot_switches(show.label = show.label,
                               add_legend = add_legend)
    norm_out <- .norm_layout_edge_type(layout, edge.type)
    layout <- norm_out$layout
    edge.type <- norm_out$edge.type
    #
    graph_data <- .get_graph_data(x)
    label_out <- .add_graph_node_labels(graph_data, show.label)
    graph_data <- label_out$df
    show.label <- label_out$show.label
    vis_out <- .incorporate_graph_vis(graph_data,
                                      se = y,
                                      edge.colour.by = edge.colour.by,
                                      edge.width.by = edge.width.by,
                                      colour_by = colour_by,
                                      shape.by = shape.by,
                                      size.by = size.by,
                                      assay.type = assay.type,
                                      other.fields = other.fields,
                                      type = type)
    graph_data <- vis_out$df
    edge.colour.by <- vis_out$edge.colour.by
    edge.width.by <- vis_out$edge.width.by
    colour_by <- vis_out$colour_by
    shape.by <- vis_out$shape.by
    size.by <- vis_out$size.by
    .graph_plotter(graph_data,
                   layout = layout,
                   edge.type = edge.type,
                   add_legend = add_legend,
                   show.label = show.label,
                   edge.colour.by = edge.colour.by,
                   edge.width.by = edge.width.by,
                   colour_by = colour_by,
                   shape.by = shape.by,
                   size.by = size.by,
                   ...)
}

################################################################################

#' @importFrom tidygraph as_tbl_graph
.get_graph_data <- function(graph){
    graph_data <- as_tbl_graph(graph)
    graph_data
}

#' @importFrom tidygraph activate
#' @importFrom dplyr mutate
.add_graph_node_labels <- function(graph_data, show.label){
    if(!("label" %in% .colnames_tbl_graph(graph_data, "nodes")) &&
       ("name" %in% .colnames_tbl_graph(graph_data, "nodes"))){
        graph_data <- graph_data %>%
            activate("nodes") %>%
            mutate(label = .data$name)
    }
    
    if(!is.logical(show.label) || length(show.label) > 1L) {
        data <- graph_data %>% 
            activate("nodes") %>%
            as_tibble()
        if(is.character(show.label) && 
                  length(show.label) == nrow(data)) {
            graph_data <- graph_data %>% 
                activate("nodes") %>%
                mutate(label = show.label)
            show.label <- TRUE
        } else if(!("label" %in% .colnames_tbl_graph(graph_data, "nodes"))){
            warning("If 'show.label' is a character vector with length != ",
                    "number of nodes in the graph or a logical/integer ",
                    "vector, a 'name' or 'label' column must exist in the ",
                    "graph data.",
                    call. = FALSE)
            show.label <- FALSE
        } else {
            if(is.numeric(show.label)){
                if(any(show.label != as.integer(show.label)) ||
                       min(show.label) < 1 ||
                       max(show.label) > nrow(data)){
                    stop("If 'show.label' is numeric, values have to be whole ",
                         "numbers and must be between 1 and the number of nodes ",
                         "in the graph",
                         call. = FALSE)
                }
                label <- rep(FALSE, nrow(data))
                label[show.label] <- TRUE
                show.label <- label
            } else if(is.character(show.label)) {
                show.label <- data$label %in% show.label
            }
            if(is.logical(show.label) &&
               length(show.label) != nrow(data)){
                stop("If 'show.label' is logical, it must have the length as ",
                     "nodes are in the graph.",
                     call. = FALSE)
            }
            graph_data <- graph_data %>% 
                activate("nodes") %>%
                mutate(label = ifelse(show.label, label, NA_character_))
            show.label <- TRUE
        }
        if(all(is.na(graph_data %>% activate("nodes") %>% pull("label")))){
            show.label <- FALSE
            warning("No labels to plot.", call. = FALSE)
        }
    }
    return(list(df = graph_data,
                show.label = show.label))
}

#' @importFrom tidygraph activate as_tibble
.colnames_tbl_graph <- function(graph_data, type){
    graph_data %>% 
        activate(!!sym(type)) %>% 
        as_tibble() %>% 
        colnames()
}

#' @importFrom tidygraph activate
.add_graph_data_or_warn <- function(data, graph_data, type, name = data$name){
    names <- .colnames_tbl_graph(graph_data, type)
    if(name %in% names){
        warning("Data for '",name,"' already present in graph nodes.",
                "Data will not be added.",
                call. = FALSE)
    }
    graph_data %>% 
        activate(!!sym(type)) %>%
        mutate(!!sym(data$name) := data$value)
}

#' @importFrom scater retrieveFeatureInfo retrieveCellInfo
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidygraph activate as_tibble
.incorporate_graph_vis <- function(graph_data,
                                   se,
                                   edge.colour.by,
                                   edge.width.by,
                                   colour_by,
                                   shape.by,
                                   size.by,
                                   assay.type = "counts",
                                   other.fields = list(),
                                   type = c("row","column")){
    type <- match.arg(type)
    type_FUN <- switch(type,
                       row = scater::retrieveFeatureInfo,
                       column = scater::retrieveCellInfo)
    variables <- c(colour_by = colour_by,
                   shape.by = shape.by,
                   size.by = size.by)
    colour_by <- NULL
    shape.by <- NULL
    size.by <- NULL
    # node data
    if(!is.null(variables)){
        # remove any variables values, which are already available and
        # rename columns by their usage
        cn <- .colnames_tbl_graph(graph_data,"nodes")
        if(length(cn) > 0L){
            f <- variables %in% cn
            if(any(f)){
                for(i in seq_along( variables[f])){
                    var_name <- names(variables[f])[i]
                    # mirror back variable name
                    assign(var_name, .get_new_var_name_value(get(var_name),
                                                             variables[f][i]))
                    # rename columns by their usage
                    graph_data %>% 
                        activate("nodes") %>%
                        dplyr::rename(!!sym(var_name) := variables[f][i])
                }
                variables <- variables[!f]
            }
        }
        if(length(variables) > 0L){
            dim_graph_nodes <- graph_data %>% 
                activate("nodes") %>%
                as_tibble() %>%
                dim()
            dim_se <- switch(type,
                             row = nrow(se),
                             column = ncol(se))
            if(dim_graph_nodes[1] != dim_se){
                stop("The number of nodes in the graph and chosen dimension ",
                     "of the SummarizedExperiment must be equal.",
                     call. = FALSE)
            }
            for(i in seq_along(variables)){
                # get data
                feature_info <- type_FUN(se, variables[i],
                                         exprs_values = assay.type)
                feature_info_name <- feature_info$name
                # mirror back variable name, if a partial match was used
                var_name <- names(variables)[i]
                assign(var_name, .get_new_var_name_value(get(var_name),
                                                         feature_info$name))
                # rename columns by their usage
                feature_info$name <- var_name
                graph_data <- .add_graph_data_or_warn(feature_info, graph_data,
                                                      type = "nodes",
                                                      feature_info_name)
            }
        }
    }
    if(length(other.fields) != 0L){
        for (o in other.fields) {
            other <- type_FUN(se, o, exprs_values = assay.type)
            graph_data <- .add_graph_data_or_warn(other, graph_data, 
                                                  type = "nodes")
        }
    }
    # edge data
    variables <- c(edge.colour.by = edge.colour.by,
                   edge.width.by = edge.width.by)
    edge.colour.by <- NULL
    edge.width.by <- NULL
    cn <- .colnames_tbl_graph(graph_data,"edges")
    variables <- variables[variables %in% cn]
    if(length(variables)  != 0L){
        for(i in seq_along(variables)){
            var_name <- names(variables)[i]
            # mirror back variable name
            assign(var_name, .get_new_var_name_value(get(var_name),
                                                     variables[i]))
            # rename columns by their usage
            graph_data <- graph_data %>% 
                activate("edges") %>%
                mutate(!!sym(var_name) := !!sym(unname(variables[i])))
        }
    }
    #
    return(list(df = graph_data,
                edge.colour.by = edge.colour.by,
                edge.width.by = edge.width.by,
                colour_by = colour_by,
                shape.by = shape.by,
                size.by = size.by))
}

.graph_plotter <- function(object,
                           layout = "kk",
                           edge.type = "fan",
                           algorithm = NULL,
                           add_legend = TRUE,
                           show.label = FALSE,
                           edge.colour.by = NULL,
                           edge.width.by = NULL,
                           colour_by = NULL,
                           shape.by = NULL,
                           size.by = NULL,
                           line_alpha = 1,
                           line_width = NULL,
                           line_width_range = c(0.5,3),
                           point.alpha = 1,
                           point.sizze = 2,
                           point_size_range = c(1,4)){
    # assemble arg list
    point_out <- .get_point_args(colour_by,
                                 shape.by,
                                 size.by,
                                 alpha = point.alpha,
                                 size = point.sizze)
    edge_out <- .get_graph_edge_args(edge.colour.by,
                                     edge.width.by,
                                     alpha = line_alpha,
                                     size = line_width,
                                     edge.type)
    edge_FUN <- match.fun(paste0("geom_edge_",edge.type))
    # begin plotting
    if(!is.null(algorithm)){
        plot_out <- ggraph(object, layout = layout, algorithm = algorithm)
    } else {
        plot_out <- ggraph(object, layout = layout)
    }
    plot_out <- plot_out +
        do.call(edge_FUN, edge_out$args) +
        do.call(geom_node_point, point_out$args)
    # add node labels
    plot_out <- .add_graph_labels(plot_out, show.label)
    # adjust edge colours
    if(!is.null(edge.colour.by)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          object %>% 
                                              activate("edges") %>% 
                                              pull("edge.colour.by"),
                                          edge.colour.by,
                                          type = "edges",
                                          na.translate = FALSE,
                                          # Specify guide
                                          guide = "edge_colourbar"
                                          )
    }
    if (!is.null(edge.width.by)) {
        if(is.numeric(object %>% activate("edges") %>% pull("edge.width.by"))){
            SIZEFUN <- scale_edge_width_continuous
        } else {
            SIZEFUN <- scale_edge_width_discrete
        }
        plot_out <- .add_extra_guide_graph(plot_out, edge.width.by) +
            SIZEFUN(range = line_width_range)
    }
    if(!is.null(size.by)){
        if(is.numeric(object %>% activate("nodes") %>% pull("size.by"))){
            SIZEFUN <- scale_size_continuous
        } else {
            SIZEFUN <- scale_size_discrete
        }
        plot_out <- plot_out +
            SIZEFUN(range = point_size_range)
    }
    # adjust point colours
    if(!is.null(colour_by)){
        plot_out <- .resolve_plot_colours(plot_out,
                                          object %>% 
                                              activate("nodes") %>% 
                                              pull("colour_by"),
                                          colour_by,
                                          fill = point_out$fill,
                                          na.translate = FALSE)
    }
    
    # add additional guides
    plot_out <- .add_extra_guide(plot_out, shape.by, size.by)
    # add theme
    plot_out <- .theme_plotGraph(plot_out)
    # optionally hide legends
    if (!add_legend) {
        plot_out <- plot_out +
            theme(legend.position = "none")
    }
    plot_out
}

#' @importFrom tidyr drop_na
.add_graph_labels <- function(plot_out, show.label){
    label <- NULL # disable note: no global binding for variable
    if(show.label){
        label_data <- plot_out$data %>% drop_na(label)
        plot_out <- plot_out +
            geom_node_label(mapping = aes(label = .data[["label"]]), 
                            data = label_data,
                            repel = TRUE,
                            max.overlaps = 100)
    }
    plot_out
}

.theme_plotGraph <- function(plot){
    plot + 
        theme_graph(base_family = "",
                    background = NA)
}
