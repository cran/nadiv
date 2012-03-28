wrap_DFC <- function(x, grandparents){
    apply(grandparents[min(x):max(x), ], MARGIN = 1, FUN = DFC)
}
