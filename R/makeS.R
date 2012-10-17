makeS <- function(pedigree, heterogametic = "0", returnS = FALSE)
{

    if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")

    damsex <- pedigree[unique(pedigree[, 2])[-1], 4]
    if(any(damsex == heterogametic)){
       pedname <- names(pedigree)
       pedigree <- pedigree[, c(1,3,2,4)]
       names(pedigree) <- pedname
       warning("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system")
    }
    numped <- numPed(pedigree[, 1:3])
    sex <- rep(-998, dim(pedigree)[1])
    sex[which(pedigree[,4] != heterogametic)] <- 1
    sex[which(pedigree[,4] == heterogametic)] <- 0
    pedin <- cbind(numped, as.integer(sex))
    N <- dim(pedin)[1]
    N2 <- N + 1


 cat(paste("Starting to make S-inverse..."))
       dnmiss <- which(pedin[,2] != -998)
       fsnmiss <- which(pedin[,3] != -998 & pedin[,4] == 1)
       Q.col <- c(pedin[,1][dnmiss], pedin[,1][fsnmiss], 1:N) 
       Q.row <- c(pedin[,2][dnmiss], pedin[,3][fsnmiss], 1:N)
       Q.x <- c(rep(-0.5, length(dnmiss)), rep(-1, length(fsnmiss)), rep(1, N))
       ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
       Q <- Matrix(0, N, N, sparse = TRUE)
       Q[1, 2] <- 1
       Q@i <- as.integer(Q.row[ord] -1)
       Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
       Q@x <- as.double(Q.x[ord])

       T <- as(solve(Q), "dgCMatrix")
       numped[numped == -998] <- N2
       Vii <- (sex + 1)/2
       f <- c(rep(0, N), -1)
       nonfound <- which(numped[,2] != N2 | numped[,3] != N2)


       Cout <- .C("vii", 
		as.integer(N), #N
		as.integer(numped[,2] - 1), #dam
		as.double(Vii), #sex
		as.integer(length(nonfound)), #nN
		as.integer(nonfound - 1), #nonfound
		as.integer(T@i), # iTP
		as.integer(c(T@p, length(T@x))), # pTP
		as.double(T@x), # xTP
		as.integer(length(T@x)), #nzmaxTP
		as.double(Vii), #v
		as.double(f)) #f



    Vinv <- Diagonal(N, 1/(Cout[[10]]))
    Sinv <- Q%*%Vinv%*%t(Q)
    listSinv <- sm2list(Sinv, rownames = as.character(pedigree[, 1]), colnames = c("Row", "Column", "Sinverse"))
    Sinv@Dimnames <- list(pedigree[,1], NULL) 
 cat(paste(".done", "\n"))
   

    if(returnS){
     cat(paste("Starting to make S..."))
       S <- as(t(T) %*% Diagonal(N, Cout[[10]]) %*% T, "dgCMatrix")
     cat(paste(".done", "\n"))

    } else S <- NULL 

return(list(S = S, Sinv = Sinv, listSinv = listSinv))
}

