rnls <- function(
formula,
data,
start,
weights = NULL,
na.action = na.fail,
psi = psi.huber,
test.vec = c("resid", "coef", "w"),
maxit = 20,
acc = 1e-06,   
trace = FALSE, 
control = nls.control(),
...)
{
    mf <- match.call()    
    formula <- as.formula(formula)
    if (length(formula) != 3)
        stop("'formula' should be a formula of the type 'y  ~ f(x, alpha)'")
    test.vec <- match.arg(test.vec)
    varNames <- all.vars(formula)
    dataName <- substitute(data)
    data <- as.data.frame(data)
    if (length(pnames <- names(start)) != length(start))
        stop("'start' must be fully named (list or numeric vector)")
    if (!((is.list(start) && all(sapply(start, is.numeric))) ||
        (is.vector(start) && is.numeric(start))) || any(is.na(match(pnames,
        varNames))))
        stop("'start' must be a list or numeric vector named with parameters in 'formula'")
    if ("w" %in% varNames || "w" %in% pnames || "w" %in% names(data))
        stop("Do not use 'w' as a variable name or as a parameter name")
    if (!is.null(weights)) {
        if (length(weights) != nrow(data))
            stop("'length(weights)' must equal the number of observations")
        if (any(weights < 0) || any(is.na(weights)))
            stop("'weights' must be nonnegative and not contain NAs")
    }
    irls.delta <- function(old, new) sqrt(sum((old - new)^2,
        na.rm = TRUE)/max(1e-20, sum(old^2, na.rm = TRUE)))
    coef <- start
    fit <- eval(formula[[3]], c(as.list(data), start))     
    y <- eval(formula[[2]], as.list(data))
    resid <- y - fit
    w <- rep(1, nrow(data))
    if (!is.null(weights))
        w <- w * weights
    oform <- formula
    formula <- as.formula(substitute(~(LHS - RHS) * w, list(LHS = formula[[2]],
        RHS = formula[[3]])))
    converged <- FALSE
    status <- "converged"
    method.exit <- FALSE
    for (iiter in 1:maxit) {
        if (trace)
            cat("robust iteration", iiter, "\n")
        previous <- get(test.vec)
        Scale <- median(abs(resid), na.rm = TRUE)/0.6745
        if (Scale == 0) {
            convi <- 0
            method.exit <- TRUE
            warning(status <- "could not compute scale of residuals")
        }
        else {
            w <- psi(resid/Scale, ...)
            if (!is.null(weights))
                w <- w * weights
            data$w <- sqrt(w)
            out <- nls(formula, data = data, start = start, algorithm = "port",
                trace = trace, na.action = na.action, control = control)
            coef <- coefficients(out)
            start <- coef
            resid <- -residuals(out)/sqrt(w)
            convi <- irls.delta(previous, get(test.vec))
        }
        converged <- convi <= acc
        if (converged)
            break
    }
    if (!converged && !method.exit)
        warning(status <- paste("failed to converge in", maxit,
            "steps"))
    if (!is.null(weights)) {
        tmp <- weights != 0
        w[tmp] <- w[tmp]/weights[tmp]
    }
    
    out <- list(m = out$m, call = match.call(), formula = oform,
        new.formula = formula, Scale = Scale, w = w,
        status = status, psi = psi, data = dataName, dataClasses = attr(attr(mf,
            "terms"), "dataClasses"))
            
    out$m$fitted <- function() fit
    out$m$lhs <- function() y
    out$call$algorithm <- "port"
    out$call$control <- control
    out$call$trace <- FALSE
    out$call$model <- TRUE
    out$call$lower <- -Inf
    out$call$upper <- Inf       
    out$message <- paste("converged in", iiter, "iterations") 
    out$control <- control       
            
    class(out) <- "nls"
    return(out)
}

