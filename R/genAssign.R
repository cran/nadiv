genAssign <- function(pedigree)
{ 
   n <- dim(pedigree)[1]
   dams <- match(pedigree[,2], pedigree[,1], nomatch = 0)
   sires <- match(pedigree[,3], pedigree[,1], nomatch = 0)
    
    if (n == 1) return(0)
    parents <- which(dams == 0 & sires == 0)
    depth <- rep(0, n)
    for (i in 1:n) {
        off <- match(dams, parents, nomatch = 0) + match(sires, 
            parents, nomatch = 0)
        if (all(off == 0)) break
        if (i == n) {
		i
               stop("Impossible pedigree: someone is their own ancestor")
	}
        parents <- which(off > 0)
        depth[parents] <- i
    }
    chaseup <- function(x, dams, sires) {
        new <- c(dams[x], sires[x])
        new <- new[new > 0]
        while (length(new) > 1) {
            x <- unique(c(x, new))
            new <- c(dams[new], sires[new])
            new <- new[new > 0]
        }
        x
    }
    dads <- sires[dams > 0 & sires > 0]
    moms <- dams[dams > 0 & sires > 0]
    dups <- duplicated(dads + moms * n)
    if (any(dups)) {
        dads <- dads[!dups]
        moms <- moms[!dups]
    }
    npair <- length(dads)
    done <- rep(FALSE, npair)
    while (TRUE) {
        pairs.to.fix <- (1:npair)[(depth[dads] != depth[moms]) & 
            !done]
        if (length(pairs.to.fix) == 0) 
            break
        temp <- pmax(depth[dads], depth[moms])[pairs.to.fix]
        who <- min(pairs.to.fix[temp == min(temp)])
        good <- moms[who]
        bad <- dads[who]
        if (depth[dads[who]] > depth[moms[who]]) {
            good <- dads[who]
            bad <- moms[who]
        }
        abad <- chaseup(bad, dams, sires)
        if (length(abad) == 1 && sum(c(dads, moms) == bad) == 
            1) {
            depth[bad] <- depth[good]
        }
        else {
            agood <- chaseup(good, dams, sires)
            tdad <- dads[-who]
            tmom <- moms[-who]
            while (1) {
                spouse <- c(tmom[!is.na(match(tdad, agood))], 
                  tdad[!is.na(match(tmom, agood))])
                temp <- unique(c(agood, spouse))
                temp <- unique(chaseup(temp, dams, sires))
                kids <- (!is.na(match(dams, temp)) | !is.na(match(sires, 
                  temp)))
                temp <- unique(c(temp, (1:n)[kids & depth <= 
                  depth[good]]))
                if (length(temp) == length(agood)) 
                  break
                else agood <- temp
            }
            if (all(match(abad, agood, nomatch = 0) == 0)) {
                depth[abad] <- depth[abad] + (depth[good] - depth[bad])
                for (i in 0:n) {
                  parents <- which(depth == i)
                  child <- match(dams, parents, nomatch = 0) + 
                    match(sires, parents, nomatch = 0)
                  if (all(child == 0)) 
                    break
                  depth[child > 0] <- pmax(i + 1, depth[child > 
                    0])
                }
            }
        }
        done[who] <- TRUE
    }
    depth
}
