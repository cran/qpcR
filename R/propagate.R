propagate <- function (fun, vals, type = c("raw", "stat"), cov = FALSE, do.sim = FALSE, cov.sim = FALSE, nsim = 10000)
{
      require(MASS, quietly = TRUE)
      type <- match.arg(type)
      expr <- as.expression(substitute(fun))
      vals <- as.matrix(vals)

      
      if (nsim > 0 && nsim < 5000) stop("Do at least 5000 simulations...")
      if (nsim < 0) stop("nsim must be >= 0!")
      
      m <- match(all.vars(expr), colnames(vals))
      if (any(is.na(m))) stop("Variable names of input dataframe and expression do not match!")
      if (length(unique(m)) != length(m)) stop("Some variable names are repetitive!")

      if (!is.logical(cov)) {
            if (!is.matrix(cov)) stop("cov must be a square covariance matrix!")
            if (is.matrix(cov)) {
                  if (dim(cov)[1] != max(m) || dim(cov)[2] != max(m)) stop(paste("cov is not a ", max(m), "x", max(m), " matrix!", sep=""))
            }
      }

      if (type == "raw") {
            meanvals <- apply(vals, 2, function(x) mean(x, na.rm = TRUE))
            sdvals <- apply(vals, 2, function(x) sd(x, na.rm = TRUE))
      } else {
            meanvals <- vals[1, ]
            sdvals <- vals[2, ]
      }
      
      if (type == "raw" && cov == TRUE) Sigma <- cov(vals, use = "complete.obs")
      if (type == "raw" && cov == FALSE) Sigma <- diag(diag(cov(vals, use = "complete.obs")))
      if (type == "stat") Sigma <- diag(vals[2, ]^2)
      
      if (is.matrix(cov)) {
            m <- match(colnames(cov), colnames(vals))
            if (any(is.na(m))) stop("Names of input dataframe and var-cov matrix do not match!")
            if (length(unique(m)) != length(m)) stop("Some names of the var-cov matrix are repetitive!")
            if (is.unsorted(m)) stop("Names of input dataframe and var-cov matrix not in the same order!")
            Sigma <- cov
      }
      
      colnames(Sigma) <- colnames(vals)
      rownames(Sigma) <- colnames(vals)
      
      if (do.sim) {
            simDat <- matrix(nrow = nsim, ncol = ncol(vals))
            
            if (!cov.sim) {
                  for (i in 1:ncol(vals)) {
                        simDat[, i] <- rnorm(nsim, mean = meanvals[i], sd = sdvals[i])
                  }
            }
            else {
                  simDat <- mvrnorm(nsim, mu = meanvals, Sigma = Sigma, empirical = TRUE)
            }
            colnames(simDat) <- colnames(vals)
      }
      else simDat <- NULL

      origDat <- t(meanvals)
      finalDat <- rbind(origDat, simDat)

      derivs <- try(lapply(colnames(vals), D, expr = expr), silent = TRUE)
      if (inherits(derivs, "try-error")) stop(paste("Error within derivs:", derivs))

      resfun <- vector(length = nrow(finalDat))
      error <- vector(length = nrow(finalDat))

      for (i in 1:nrow(finalDat)) {
            nderivs <- sapply(derivs, eval, envir = as.list(finalDat[i, ]))
            resfun[i] <- eval(expr, envir = as.list(finalDat[i, ]))
            error[i] <- c(nderivs %*% Sigma %*% matrix(nderivs))
      }
      if(do.sim && length(unique(error)) == 1) print("Monte Carlo simulation gave unique repetitive error value! Are all derivations constants?")

      return(list(errProp = sqrt(error[1]), evalProp = resfun[1], errSim = ifelse(do.sim, mean(sqrt(error[-1])), NA),
                  evalSim = ifelse(do.sim, mean(resfun[-1]), NA), derivs = derivs, covMat = Sigma, simVec = error[-1]))
}