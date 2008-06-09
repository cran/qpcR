ratiocalc.pcrbatch <- function (data, group = NULL, ratio = c("ind", "first"), which.eff = c("sig", "sli", "exp"),
                                    iter = c("combs", "perms"), rep.all = TRUE, ttest = c("cp", "Ecp"), ...)
{
      require(gtools, quietly = TRUE)
    
      if (class(data) != "pcrbatch")
            stop("data is not of class 'pcrbatch'!")
      if (!is.numeric(ratio))  {
            ratio <- match.arg(ratio)
            RATIO <- NULL
      }
      else {
            if (ratio < 1 || ratio > 2) stop("Efficiency should be between 1 and 2!")
            RATIO <- ratio
            ratio <- "NUMERIC"
      }
      
      which.eff <- match.arg(which.eff)
      iter <- match.arg(iter)
      ttest <- match.arg(ttest)

      if (!is.null(group)) {
            if (is.factor(group)) group <- as.numeric(group)
            if (length(group) != ncol(data) - 1) stop("Length of 'group' and 'data' do not match!")
            exp <- which(group < 100)
            con <- which(group >= 100)
            EXP <- as.factor(group[exp])
            CON <- as.factor(group[con])
            if (length(con) > 0) {
                  CONYES <- TRUE
                  if (length(exp) != length(con))
                        stop("Number of target PCRs and reference PCRs do not match!")
                  if (nlevels(EXP) != nlevels(CON))
                        stop("Some target PCRs have no reference PCR!")
                  if (nlevels(EXP) == 1)
                        stop("Ratio calculation not possible with one target PCR!")
            }
            else CONYES <- FALSE
      }
      else {
            exp <- 1:(ncol(data) - 1)
            CONYES <- FALSE
            EXP <- as.factor(exp)
            group <- as.factor(exp)
      }
    
      rn <- data[, 1]
      sel.eff <- switch(which.eff, sig = "sig.eff", sli = "sli.effmax",
                        exp = "exp.eff")
      ALL.eff <- as.numeric(data[which(rn == sel.eff), -1])
      ALL.cp <- as.numeric(data[which(rn == "sig.cpD2"), -1])
      EXP.eff <- ALL.eff[exp]
      EXP.eff <- split(EXP.eff, EXP)
      EXP.cp <- ALL.cp[exp]
      EXP.cp <- split(EXP.cp, EXP)
    
      if (CONYES) {
            CON.eff <- ALL.eff[con]
            CON.eff <- split(CON.eff, CON)
            CON.cp <- ALL.cp[con]
            CON.cp <- split(CON.cp, CON)
      }
    
      len1 <- sapply(EXP.eff, function(x) length(x))
      len2 <- sapply(EXP.cp, function(x) length(x))
    
      if (CONYES) {
        len3 <- sapply(CON.eff, function(x) length(x))
        len4 <- sapply(CON.cp, function(x) length(x))
      }
    
      lenAll <- c(len1, len2, if (exists("len3")) len3, if (exists("len4")) len4)
      mlen <- max(lenAll)
      expEff <- as.data.frame(lapply(EXP.eff, function(x) x <- c(x, rep(NA, mlen - length(x)))))
      expCp <- as.data.frame(lapply(EXP.cp, function(x) x <- c(x, rep(NA, mlen - length(x)))))
    
      if (CONYES) {
            conEff <- as.data.frame(lapply(CON.eff, function(x) x <- c(x, rep(NA, mlen - length(x)))))
            conCp <- as.data.frame(lapply(CON.cp, function(x) x <- c(x, rep(NA, mlen - length(x)))))
      }
      else {
            conEff <- NA
            conCp <- NA
      }
    
      effDat <- cbind(expEff, conEff)
      cpDat <- cbind(expCp, conCp)
    
      if (iter == "combs") COMBS <- combinations(ncol(expEff), 2, repeats.allowed = rep.all)
            else COMBS <- permutations(ncol(expEff), 2, repeats.allowed = rep.all)
            
      OUT <- NULL
      listNAMES <- colnames(data)[-1]
      groupNAMES <- split(listNAMES, as.factor(group))
      allNAMES <- NULL
      PROPLIST <- list()
      n <- NULL
    
      EXPRESSIONS <- list(
                        EXPR1 = expression((E2^cp2)/(E1^cp1)),
                        EXPR2 = expression(((E2^cp2)/(E1^cp1))/((E4^cp4)/(E3^cp3))),
                        EXPR3 = expression((E1^cp2)/(E1^cp1)),
                        EXPR4 = expression(((E1^cp2)/(E1^cp1))/((E3^cp4)/(E3^cp3))),
                        EXPR5 = as.expression(as.list(substitute(expression((val^cp2)/(val^cp1)), list(val = RATIO)))[-1]),
                        EXPR6 = as.expression(as.list(substitute(expression(((val^cp2)/(val^cp1))/(val^cp4)/(val^cp3)), list(val = RATIO)))[-1])
      )

      HTESTS <- list(
                  HTEST1 = expression(t.test(cp1, cp2, ...)),
                  HTEST2 = expression(t.test(cp1 - cp3, cp2 - cp4, ...)),
                  HTEST3 = expression(t.test(E1^cp1, E2^cp2, ...)),
                  HTEST4 = expression(t.test((E1^cp1) - (E2^cp2), (E3^cp3) - (E4^cp4), ...)),
                  HTEST5 = expression(t.test(E1^cp1, E1^cp2, ...)),
                  HTEST6 = expression(t.test((E1^cp1) - (E1^cp2), (E3^cp3) - (E3^cp4), ...)),
                  HTEST7 = as.expression(as.list(substitute(expression(t.test(val^cp1, val^cp2, ...)), list(val = RATIO)))[-1]),
                  HTEST8 = as.expression(as.list(substitute(expression(t.test((val^cp1) - (val^cp2), (val^cp3) - (val^cp4), ...)) ,list(val = RATIO)))[-1])
       )

      CHOICES <- expand.grid(ratio = c("ind", "first", "NUMERIC"), CONYES = c(FALSE, TRUE), ttest = c("cp", "Ecp"))
      expr <- c(1, 3, 5, 2, 4, 6, 1, 3, 5, 2, 4, 6)
      htest <- c(1, 1, 1, 2, 2, 2, 3, 5, 7, 4, 6, 8)
      CHOICE <- c(ratio, CONYES, ttest)
      WHICH <- sapply(CHOICE, function(x) which(x == CHOICES, arr.ind = TRUE)[, 1])
      TABLE <- as.data.frame(table(unlist(WHICH)))
      IND <- as.numeric(as.vector(TABLE$Var1[TABLE$Freq == length(WHICH)]))

      EXPR <- EXPRESSIONS[[expr[IND]]]
      HTEST <- HTESTS[[htest[IND]]]

      for (i in 1:nrow(COMBS)) {
            whichcol <- COMBS[i, ]
            E1 <- expEff[, whichcol[1]]
            E2 <- expEff[, whichcol[2]]
            cp1 <- expCp[, whichcol[1]]
            cp2 <- expCp[, whichcol[2]]
            expFrame <- cbind(E1, E2, cp1, cp2)
            expconFrame <- expFrame

            if (CONYES) {
                  E3 <- conEff[, whichcol[1]]
                  E4 <- conEff[, whichcol[2]]
                  cp3 <- conCp[, whichcol[1]]
                  cp4 <- conCp[, whichcol[2]]
                  conFrame <- cbind(E3, E4, cp3, cp4)
                  expconFrame <- cbind(expFrame, conFrame)
            }

            PROP <- propagate(EXPR, expconFrame, ...)
            STAT <- try(eval(HTEST), silent = TRUE)
            if (inherits(STAT, "try-error")) STAT <- list(p.value = -1)

            OUTtemp <- unlist(PROP[1:6])
            nobs <- nrow(expconFrame)
            OUTtemp <- c(OUTtemp, n = nobs, t.test = STAT$p.value)
            OUT <- cbind(OUT, OUTtemp)
            firstNAME <- groupNAMES[[COMBS[i, 1]]][1]
            secondNAME <- groupNAMES[[COMBS[i, 2]]][1]
            NAME <- paste(firstNAME, "/", secondNAME, sep = "")
            allNAMES <- c(allNAMES, NAME)
            PROPLIST[[i]] <- PROP
      }
    
      colnames(OUT) <- allNAMES
      OUTclip <- cbind(rownames(OUT), OUT)
      write.table(OUTclip, file = "clipboard-64000", sep = "\t",
                  row.names = FALSE)
      RET <- list(ratios = OUT, propList = PROPLIST)
      class(RET) <- "ratiocalc"
      invisible(RET)
}


