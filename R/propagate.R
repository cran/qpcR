propagate <- function (expr, data, type = c("raw", "stat"), cov = FALSE, do.sim = FALSE, cov.sim = FALSE, nsim = 10000)
{
      require(MASS, quietly = TRUE)
      type <- match.arg(type)
      if (!is.expression(expr) && !is.call(expr)) stop("'expr' must be an expression")
      DATA <- as.matrix(data)
      EXPR <- expr

      
      if (nsim > 0 && nsim < 5000) stop("Do at least 5000 simulations...")
      if (nsim < 0) stop("nsim must be >= 0!")
      
      m <- match(all.vars(expr), colnames(data))
      if (any(is.na(m))) stop("Variable names of input dataframe and expression do not match!")
      if (length(unique(m)) != length(m)) stop("Some variable names are repetitive!")

      if (!is.logical(cov)) {
            if (!is.matrix(cov)) stop("cov must be a square covariance matrix!")
            if (is.matrix(cov)) {
                  if (dim(cov)[1] != max(m) || dim(cov)[2] != max(m)) stop(paste("cov is not a ", max(m), "x", max(m), " matrix!", sep=""))
            }
      }

      if (type == "raw") {
            meanvals <- apply(DATA, 2, function(x) mean(x, na.rm = TRUE))
            sdvals <- apply(DATA, 2, function(x) sd(x, na.rm = TRUE))
      } else {
            meanvals <- DATA[1, ]
            sdvals <- DATA[2, ]
      }
      
      if (type == "raw" && cov == TRUE) SIGMA <- cov(DATA, use = "complete.obs")
      if (type == "raw" && cov == FALSE) SIGMA <- diag(diag(cov(DATA, use = "complete.obs")))
      if (type == "stat") SIGMA <- diag(DATA[2, ]^2)
      
      if (is.matrix(cov)) {
            m <- match(colnames(cov), colnames(DATA))
            if (any(is.na(m))) stop("Names of input dataframe and var-cov matrix do not match!")
            if (length(unique(m)) != length(m)) stop("Some names of the var-cov matrix are repetitive!")
            if (is.unsorted(m)) stop("Names of input dataframe and var-cov matrix not in the same order!")
            SIGMA <- cov
      }
      
      colnames(SIGMA) <- colnames(DATA)
      rownames(SIGMA) <- colnames(DATA)
      
      if (do.sim) {
            simDat <- matrix(nrow = nsim, ncol = ncol(DATA))
            
            if (!cov.sim) {
                  for (i in 1:ncol(DATA)) {
                        simDat[, i] <- rnorm(nsim, mean = meanvals[i], sd = sdvals[i])
                  }
            }
            else {
                  simDat <- mvrnorm(nsim, mu = meanvals, Sigma = SIGMA, empirical = TRUE)
            }
            colnames(simDat) <- colnames(DATA)
      }
      else simDat <- NULL

      origDat <- t(meanvals)
      finalDat <- rbind(origDat, simDat)
      
      if (type == "raw") {
            rawDat <- data[complete.cases(data), ]
            if (!is.matrix(rawDat)) rawDat <- t(rawDat)
      }

      derivs <- try(lapply(colnames(DATA), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))

      resfun <- NULL
      error <- NULL

      for (i in 1:nrow(finalDat)) {
            nderivs <- sapply(derivs, eval, envir = as.list(finalDat[i, ]))
            resfun[i] <- eval(EXPR, envir = as.list(finalDat[i, ]))
            error[i] <- c(nderivs %*% SIGMA %*% matrix(nderivs))
      }

      if (type == "raw") {
            rawfun <- NULL
            for (i in 1:nrow(rawDat)) {
                  rawfun[i] <- eval(EXPR, envir = as.list(rawDat[i, ]))
            }
      } else rawfun <- 0
      
      if(do.sim && length(unique(error)) == 1) print("Monte Carlo simulation gave unique repetitive error value! Are all derivations constants?")

      return(list(evalExpr = resfun[1], evalSim = ifelse(do.sim, mean(resfun[-1], na.rm = TRUE), NA), errProp = sqrt(error[1]),
                  errPropSim = ifelse(do.sim, mean(sqrt(error[-1]), na.rm = TRUE), NA), errEval = sd(rawfun, na.rm = TRUE),
                  errEvalSim = ifelse(do.sim, sd(resfun[-1], na.rm = TRUE), NA), derivs = derivs, covMat = SIGMA,
                  errVec = error[-1], evalVec = resfun[-1]))
}