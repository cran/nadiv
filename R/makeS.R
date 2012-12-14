makeS <- function(pedigree, heterogametic = "0", DosageComp = c(NULL, "ndc", "hoin"), returnS = FALSE){

    if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")
    numped <- numPed(pedigree[, 1:3])

    damsex <- pedigree[unique(numped[, 2])[-1], 4]
    if(any(damsex == heterogametic)){
       pedname <- names(pedigree)
       pedigree <- pedigree[, c(1,3,2,4)]
       names(pedigree) <- pedname
       numped <- numPed(pedigree[, 1:3])
      warning("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system")
    }
    sex <- rep(-998, dim(pedigree)[1])
    sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
    sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
    pedin <- cbind(numped, as.integer(sex))
    N <- dim(pedin)[1]
    N2 <- N + 1
    dc.model <- match.arg(DosageComp)
    if(is.null(dc.model)) dc.model <- "ndc"

    if(dc.model == "ndc"){

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
                as.double(0.25), #DC 
		as.integer(N), #N
		as.integer(numped[, 2] - 1), #dam
                as.integer(numped[, 3] - 1), #sire
		as.double(Vii), #sex
		as.integer(length(nonfound)), #nN
		as.integer(nonfound - 1), #nonfound
		as.integer(T@i), # iTP
		as.integer(c(T@p, length(T@x))), # pTP
		as.double(T@x), # xTP
		as.integer(length(T@x)), #nzmaxTP
		as.double(Vii), #v
		as.double(f)) #f

    }

    
    if(dc.model == "hoin"){
      cat(paste("Starting to make S-inverse..."))
          fdnmiss <- which(pedin[,2] != -998 & pedin[,4] == 1)
          mdnmiss <- which(pedin[,2] != -998 & pedin[,4] == 0)
          fsnmiss <- which(pedin[,3] != -998 & pedin[,4] == 1)
          Q.col <- c(pedin[,1][fdnmiss], pedin[,1][mdnmiss], pedin[,1][fsnmiss], 1:N) 
          Q.row <- c(pedin[,2][fdnmiss], pedin[,2][mdnmiss], pedin[,3][fsnmiss], 1:N)
          Q.x <- c(rep(-0.5, length(fdnmiss)), rep(-1, length(mdnmiss)), rep(-0.5, length(fsnmiss)), rep(1, N))
          ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
          Q <- Matrix(0, N, N, sparse = TRUE)
          Q[1, 2] <- 1
          Q@i <- as.integer(Q.row[ord] -1)
          Q@p <- as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1)
          Q@x <- as.double(Q.x[ord])

          T <- as(solve(Q), "dgCMatrix")
          numped[numped == -998] <- N2
          Vii <- (2 - sex)
          f <- c(rep(0, N), -1)
          nonfound <- which(numped[,2] != N2 | numped[,3] != N2)


          Cout <- .C("vii",
                as.double(1), #DC 
		as.integer(N), #N
		as.integer(numped[,2] - 1), #dam
                as.integer(numped[, 3] - 1), #sire
		as.double(Vii), #sex
		as.integer(length(nonfound)), #nN
		as.integer(nonfound - 1), #nonfound
		as.integer(T@i), # iTP
		as.integer(c(T@p, length(T@x))), # pTP
		as.double(T@x), # xTP
		as.integer(length(T@x)), #nzmaxTP
		as.double(Vii), #v
		as.double(f)) #f


    }


    Vinv <- Diagonal(N, 1/(Cout[[12]]))
    Sinv <- Q %*% Vinv %*% t(Q)
    listSinv <- sm2list(Sinv, rownames = as.character(pedigree[, 1]), colnames = c("Row", "Column", "Sinverse"))
    Sinv@Dimnames <- list(pedigree[,1], NULL) 
    cat(paste(".done", "\n"))

    if(returnS){
       cat(paste("Starting to make S..."))
          S <- as(t(T) %*% Diagonal(N, Cout[[12]]) %*% T, "dgCMatrix")
       cat(paste(".done", "\n"))
    } else{
         S <- NULL
      }

return(list(model = dc.model, S = S, Sinv = Sinv, listSinv = listSinv))
}



