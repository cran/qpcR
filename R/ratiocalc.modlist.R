ratiocalc.modlist <- function(data, group = NULL, ratio.fun = c("E1-E2", "E1-E1", "2-2"), ...)
{
      require(gtools, quietly = TRUE)
      if (class(data) != "modlist") stop("data is not of class 'modlist'!")
      ratio.fun <- match.arg(ratio.fun)

      if (!is.null(group)) {
            if (is.factor(group)) group <- as.numeric(group)
            if (length(group) != length(data)) stop("Length of 'group' and 'data' do not match!")
            exp <- which(group < 100)
            con <- which(group >= 100)
            EXP <- as.factor(group[exp])
            CON <- as.factor(group[con])
            if (length(con) > 0) {
                  CONYES <- TRUE
                  if (length(exp) != length(con)) stop("Number of target PCRs and reference PCRs do not match!")
                  if (nlevels(EXP) != nlevels(CON)) stop("Some target PCRs have no reference PCR!")
                  if (nlevels(EXP) == 1) stop("Ratio calculation not possible with one target PCR!")                  
            } else CONYES <- FALSE
      } else {
            exp <- 1:length(data)
            CONYES <- FALSE
            EXP <- as.factor(exp)
      }

      EXP.eff <- sapply(data[exp], function(x) efficiency(x, plot = FALSE)$eff)
      EXP.eff <- split(EXP.eff, EXP)
      EXP.cp <- sapply(data[exp], function(x) efficiency(x, plot = FALSE)$cpD2)
      EXP.cp <- split(EXP.cp, EXP)

      if (CONYES) {
            CON.eff <- sapply(data[con], function(x) efficiency(x, plot = FALSE)$eff)
            CON.eff <- split(CON.eff, CON)
            CON.cp <- sapply(data[con], function(x) efficiency(x, plot = FALSE)$cpD2)
            CON.cp <- split(CON.cp, CON)
      }

      len1 <- sapply(EXP.eff, function(x) length(x))
      len2 <- sapply(EXP.cp, function(x) length(x))
      
      if (CONYES) {
            len3 <- sapply(CON.eff, function(x) length(x))
            len4 <- sapply(CON.cp, function(x) length(x))
      }
      
      lenAll <- c(len1, len2, if (exists("len3")) len3, if(exists("len4")) len4)
      mlen <- max(lenAll)

      expEff <- as.data.frame(lapply(EXP.eff, function(x) x <- c(x, rep(NA, mlen - length(x)))))
      expCp <- as.data.frame(lapply(EXP.cp, function(x) x <- c(x, rep(NA, mlen - length(x)))))


      if (CONYES) {
            conEff <- as.data.frame(lapply(CON.eff, function(x) x <- c(x, rep(NA, mlen - length(x)))))
            conCp <- as.data.frame(lapply(CON.cp, function(x) x <- c(x, rep(NA, mlen - length(x)))))
      } else {
            conEff <- NA
            conCp <- NA
      }


      effDat <- cbind(expEff, conEff)
      cpDat <- cbind(expCp, conCp)
      
      COMBS <- combinations(ncol(expEff), 2)

      OUT <- NULL
      listNAMES <- sapply(data, function(x) x$names)
      allNAMES <- NULL
      PROPLIST <- list()
      
      for (i in 1:nrow(COMBS)) {
            whichcol <- COMBS[i, ]
            E1 <- expEff[, whichcol[1]]
            E2 <- expEff[, whichcol[2]]
            cp1 <- expCp[, whichcol[1]]
            cp2 <- expCp[, whichcol[2]]
            expFrame <- cbind(E1, E2, cp1, cp2)
            expconFrame <- expFrame
            
            if (CONYES)  {
                  E3 <- conEff[, whichcol[1]]
                  E4 <- conEff[, whichcol[2]]
                  cp3 <- conCp[, whichcol[1]]
                  cp4 <- conCp[, whichcol[2]]
                  conFrame <- cbind(E3, E4, cp3, cp4)
                  expconFrame <- cbind(expFrame, conFrame)
            }

            if (ratio.fun == "E1-E2") {
                  if (!CONYES) PROP <- propagate((E2^cp2)/(E1^cp1), vals = expconFrame, type = "raw", ...)
                  if (CONYES) PROP <- propagate(((E2^cp2)/(E1^cp1))/((E4^cp4)/(E3^cp3)), vals = expconFrame, type = "raw", ...)
            }
            if (ratio.fun == "E1-E1") {
                  if (!CONYES) PROP <- propagate((E1^cp1)/(E1^cp2), vals = expconFrame, type = "raw", ...)
                  if (CONYES) PROP <- propagate(((E1^cp1)/(E1^cp2))/((E3^cp3)/(E3^cp4)), vals = expconFrame, type = "raw", ...)
            }
            if (ratio.fun == "2-2") {
                  if (!CONYES) PROP <- propagate((2^cp1)/(2^cp2), vals = expconFrame, type = "raw", ...)
                  if (CONYES) PROP <- propagate(((2^cp1)/(2^cp2))/((2^cp3)/(2^cp4)), vals = expconFrame, type = "raw", ...)
            }
            
            OUTtemp <- unlist(PROP[c(2, 1, 4, 3)])
            OUT <- cbind(OUT, OUTtemp)
            NAME <-  paste(listNAMES[COMBS[i, 1]], "/", listNAMES[COMBS[i, 2]], sep = "")
            allNAMES <- c(allNAMES, NAME)
            PROPLIST[[i]] <- PROP
      }
      colnames(OUT) <- allNAMES
      OUTclip <- cbind(rownames(OUT), OUT)
      write.table(OUTclip, file = "clipboard-64000", sep = "\t", row.names = FALSE)
      invisible(list(ratios = OUT, propList = PROPLIST))
}