vchoose <- function(r, N, which2, values2){
   vapply(1:N, FUN = choose, FUN.VALUE = vector("numeric", 1), which = which2[r, ], values = values2)
}

