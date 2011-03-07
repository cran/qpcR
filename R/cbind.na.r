cbind.na <- function (..., deparse.level = 1)
{
    na <- nargs() - (!missing(deparse.level))
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)
    argl <- list(...)
      
    while (na > 0 && is.null(argl[[na]])) {
        argl <- argl[-na]
        na <- na - 1
    }
    
    ### determine max length/nrow
    tempLEN <- integer(length = length(argl))
    for (i in 1:length(argl)) {
      if (is.null(dim(argl[[i]]))) tempLEN[i] <- length(argl[[i]]) else tempLEN[i] <- nrow(argl[[i]])
    }    
    maxLEN <- max(tempLEN, na.rm = TRUE)     
    ### added NA fill to max length/nrow
    for (i in 1:length(argl)) {
      if (is.null(dim(argl[[i]]))) argl[[i]] <- c(argl[[i]], rep(NA, maxLEN - length(argl[[i]])))  
      else argl[[i]] <- rbind(argl[[i]], matrix(nrow = maxLEN - nrow(argl[[i]]), ncol = ncol(argl[[i]])))
    }    
    
    if (na == 0)
        return(NULL)
    if (na == 1) {
        if (isS4(..1))
            return(cbind2(..1))
        else return(.Internal(cbind(deparse.level, ...)))
    }
    if (deparse.level) {
        symarg <- as.list(sys.call()[-1L])[1L:na]
        Nms <- function(i) {
            if (is.null(r <- names(symarg[i])) || r == "") {
                if (is.symbol(r <- symarg[[i]]) || deparse.level ==
                  2)
                  deparse(r)
            }
            else r
        }
    }
    if (na == 2) {
        ### changed to second argl item
        r <- argl[[2]]
        #r <- ..2
        fix.na <- FALSE
    }
    else {
        nrs <- unname(lapply(argl, nrow))
        iV <- sapply(nrs, is.null)
        fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
        if (fix.na) {
            nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
            argl[[na]] <- cbind(rep(argl[[na]], length.out = nr),
                deparse.level = 0)
        }
        if (deparse.level) {
            if (fix.na)
                fix.na <- !is.null(Nna <- Nms(na))
            if (!is.null(nmi <- names(argl)))
                iV <- iV & (nmi == "")
            ii <- if (fix.na)
                2:(na - 1)
            else 2:na
            if (any(iV[ii])) {
                for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
                  names(argl)[i] <- nmi
            }
        }
        r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
    }
    d2 <- dim(r)
    ### changed to first argl item
    r <- cbind2(argl[[1]], r)
    #r <- cbind2(..1, r)
    if (deparse.level == 0)
        return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
    if (ism1 && ism2)
        return(r)
    Ncol <- function(x) {
        d <- dim(x)
        if (length(d) == 2L)
            d[2L]
        else as.integer(length(x) > 0L)
    }
    nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
    nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
    if (nn1 || nn2 || fix.na) {
        if (is.null(colnames(r)))
            colnames(r) <- rep.int("", ncol(r))
        setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams))
            ""
        else nams
        if (nn1)
            setN(1, N1)
        if (nn2)
            setN(1 + l1, N2)
        if (fix.na)
            setN(ncol(r), Nna)
    }
    r
}


