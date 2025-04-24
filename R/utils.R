
################################################################################
# internal methods loaded from other packages

.check_assay_present <- mia:::.check_assay_present
.require_package <- mia:::.require_package
.check_taxonomic_rank <- mia:::.check_taxonomic_rank
.check_for_taxonomic_data_order <- mia:::.check_for_taxonomic_data_order
.calc_rel_abund <- mia:::.calc_rel_abund

.is_a_bool <- mia:::.is_a_bool
.is_non_empty_character <- mia:::.is_non_empty_character
.is_non_empty_string <- mia:::.is_non_empty_string
.is_a_string <- mia:::.is_a_string
.are_whole_numbers <- mia:::.are_whole_numbers
.is_numeric_string <- mia:::.is_numeric_string
.is_function <- mia:::.is_function
.get_name_in_parent <- mia:::.get_name_in_parent
.is_an_integer <- mia:::.is_an_integer
TAXONOMY_RANKS <- mia:::TAXONOMY_RANKS
.is_a_numeric <- mia:::.is_a_numeric
.capitalize <- mia:::.capitalize

.norm_label <- function(label, x){
    if(!is.null(label)){
        if(is.numeric(label)){
            n_v <- seq_len(nrow(x))
            if(!all(label %in% n_v)){
                stop("If 'label' is numeric, all values must be between 1 ",
                    "and nrow(x). If rank is not NULL, the dimension might ",
                    "change.", call. = FALSE)
            }
            label <- n_v %in% label
        } else if(is.character(label)){
            if(!all(label %in% rownames(x))){
                stop("If 'label' is character, all values must be in ",
                    "rownames(x). If rank is not NULL, the rownames might ",
                    "change.", call. = FALSE)
            }
            label <- rownames(x) %in% label
        } else if(is.logical(label)){
            if(length(label) != nrow(x)){
                stop("If 'label' is logical, length(label) == nrow(x) mut be ",
                    "TRUE. If rank is not NULL, the rownames might ",
                    "change.", call. = FALSE)
            }
        } else {
            stop("'label' must be a vector.", call. = FALSE)
        }
    }
    label
}

# This function checks whether variable can be found from colData or rowData.
.check_metadata_variable <- function(
        tse, var, row = FALSE, col = FALSE, multiple = FALSE,
        var.name = .get_name_in_parent(var)){
    if( !.is_a_bool(multiple) ){
        stop("'multiple' must be TRUE or FALSE.", call. = FALSE)
    }
    # If the variable is not NULL
    if( !is.null(var) ){
        # It must be a string and found from colData/rowData
        is_string <- ifelse(multiple, is.character(var), .is_a_string(var))
        check_values <- c()
        check_values <- c(check_values, if(col) colnames(colData(tse)))
        check_values <- c(check_values, if(row) colnames(rowData(tse)))
        var_found <- all( var %in% check_values )
        if( !(is_string && var_found) ){
            stop("'", var.name, "' must be", ifelse(multiple, "", "a single "),
                "character value from the following options: '",
                paste0(check_values, collapse = "', '"), "'", call. = FALSE)
        }
    }
    return(NULL)
}
