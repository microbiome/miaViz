




#' @rdname plotAbundance
setGeneric("plotComposition", signature = c("x"),
           function(x, ...)
               standardGeneric("plotComposition"))







# .get_first_numeric_or_factor <- function(x){
#     f <- vapply(colData(x),is.numeric,logical(1)) |
#         vapply(colData(x),is.factor,logical(1))
#     if(!any(f)){
#         stop("No numeric of factor values found in colData(x).", call. = FALSE)
#     }
#     colnames(colData(x))[f][1L]
# }
# 
# .get_first_factor_or_character <- function(x){
#     f <- vapply(colData(x),is.character,logical(1)) |
#         vapply(colData(x),is.factor,logical(1))
#     if(!any(f)){
#         stop("No character of factor values found in colData(x).",
#              call. = FALSE)
#     }
#     colnames(colData(x))[f][1L]
# }
# 
# .norm_across <- function(across, x){
#     if(!is.null(across)){
#         if(!is.character(across) || length(across) == 0L){
#             stop("'across' must be a single character value.", call. = FALSE)
#         }
#         if(!(across %in% colnames(colData(x)))){
#             stop("'across' must define a column name of colData(x)",
#                  call. = FALSE)
#         }
#         values <- colData(x)[,across]
#         if(!is.numeric(values) && !is.factor(values)){
#             stop("values defined by 'across' must be numeric or a factor",
#                  call. = FALSE)
#         }
#     } else {
#         across <- .get_first_numeric_or_factor(x)
#     }
#     across
# }
# 
# .norm_per <- function(per, x){
#     if(!is.null(per)){
#         if(!is.character(per) || length(per) == 0L){
#             stop("'per' must be a single character value.", call. = FALSE)
#         }
#         if(!(per %in% colnames(colData(x)))){
#             stop("'per' must define a column name of colData(x)",
#                  call. = FALSE)
#         }
#         values <- colData(x)[,per]
#         if(!is.character(values) && !is.factor(values)){
#             stop("values defined by 'across' must be character or a factor",
#                  call. = FALSE)
#         }
#     } else {
#         per <- .get_first_factor_or_character(x)
#     }
#     per
# }