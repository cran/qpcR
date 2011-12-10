rm(list=ls(all=TRUE))
#####################
l7 <- list(
      expr = "Fluo ~ c + (k1 * Cycles) + (k2 * Cycles^2) + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)",
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            c + (k1 * x) + (k2 * x^2) + (d - c)/((1 + exp(b * (log(x) - log(e))))^f)
      },
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ log(x2))
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- exp(-coefVec[1]/b)
            f <- 1
            lmFit2 <- rlm(y2[1:10] ~ x2[1:10])
            k1 <- coef(lmFit2)[2]
            k2 <- - 0.01 * k1
            ssVal <- as.numeric(c(b, c, d, e, f, k1, k2))
            names(ssVal) <- l7$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            k1 + 2 * k2 * x + b * (c - d) * e^-b * f * x^(-1 + b) * (1 + e^-b * x^b)^(-1 - f)
      },
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            (e^(-2 * b) * (1+e^-b * x^b)^(-2 - f) * (-b * (c - d) * f * x^b * (e^b + x^b) + 2 * k2 * x^2 *
            (e^b + x^b)^2 * (1 + e^-b * x^b)^f + b^2 * (c - d) * f * x^b * (e^b - f * x^b)))/x^2
      },
      inv = function(y, parm) {
            x <- 1:100
            fn <- function(x, parm) l7$fct(x, parm) - y 
            uniroot(fn, interval = c(1, 100), parm)$root
      },
      expr.grad = expression(c + (k1 * Cycles) + (k2 * Cycles^2) + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)),
      inv.grad =  NULL,
      parnames = c("b", "c", "d", "e", "f", "k1", "k2"),
      name = "l7",
      type = "seven-parameter log-logistic"
)

l6 <- list(
      expr = "Fluo ~ c + (k * Cycles) + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)",
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            c + (k * x) + (d - c)/((1 + exp(b * (log(x) - log(e))))^f)
      },
      ssFct = function (x, y) {
            d <- max(y) + 0.001                                     
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ log(x2))
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- exp(-coefVec[1]/b)
            f <- 1
            lmFit2 <- rlm(y2[1:10] ~ x2[1:10])
            k <- coef(lmFit2)[2]
            ssVal <- as.numeric(c(b, c, d, e, f, k))
            names(ssVal) <- l6$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            k + b * (c - d) * e^-b * f * x^(-1 + b) * (1 + e^-b * x^b)^(-1-f)      
      },
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            -b * (c - d) * e^-(2 * b) * f * x^(-2 + b) * (1 + e^-b * x^b)^(-2-f) * (-(-1 + b) * e^b + (1 + b * f) * x^b)   
      },
      inv = function(y, parm) {
            x <- 1:100
            fn <- function(x, parm) l6$fct(x, parm) - y 
            uniroot(fn, interval = c(1, 100), parm)$root    
      },
      expr.grad = expression(c + (k * Cycles) + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)),
      inv.grad =  NULL,
      parnames = c("b", "c", "d", "e", "f", "k"),
      name = "l6",
      type = "six-parameter log-logistic"
)

l5 <- list(
      expr = "Fluo ~ c + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)", 
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            c + (d - c)/((1 + exp(b * (log(x) - log(e))))^f)
      },  
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))            
            lmFit <- rlm(logitTrans ~ log(x2))
            coefVec <- coef(lmFit)              
            b <- coefVec[2]
            e <- exp(-coefVec[1]/b)
            f <- 1
            ssVal <- as.numeric(c(b, c, d, e, f))
            names(ssVal) <- l5$parnames
            return(ssVal)
      },     
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            b * (c - d) * e^-b * f * x^(-1 + b) * (1 + e^-b * x^b)^(-1-f)
      },        
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            -b * (c - d) * e^(-2 * b) * f * x^(-2 + b) *
            (1 + e^-b * x^b)^(-2 - f) * (-(-1 + b) *  e^b + (1 + b * f) * x^b) 
      },  
      inv = function(y, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            e * (1/(-1 + ((c - d)/(c - y))^(1/f)))^(-1/b)
      },
      expr.grad = expression(c + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)),
      inv.grad =  expression(e * (1/(-1 + ((c - d)/(c - Fluo))^(1/f)))^(-1/b)),      
      parnames = c("b", "c", "d", "e", "f"),
      name = "l5",
      type = "five-parameter log-logistic"
)

l4 <- list(
      expr = "Fluo ~ c + (d - c)/(1 + exp(b * (log(Cycles) - log(e))))", 
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            c + (d - c)/(1 + exp(b * (log(x) - log(e))))
      },   
      ssFct = function (x, y) {            
            d <- max(y, na.rm = TRUE) + 0.01
            c <- min(y, na.rm = TRUE) - 0.01
            x2 <- x[y > 0]            
            y2 <- y[y > 0]             
            logitTrans <- log((d - y2)/(y2 - c))            
            lmFit <- rlm(logitTrans ~ log(x2))
            coefVec <- coef(lmFit)               
            b <- coefVec[2]             
            e <- exp(-coefVec[1]/b)            
            ssVal <- as.numeric(c(b, c, d, e))
            names(ssVal) <- l4$parnames
            return(ssVal)
      },   
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            (b * (c - d) * e^b * x^(-1 + b))/(e^b + x^b)^2
      },      
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            (b * (c - d) * e^b * x^(-2 + b) * ((-1 + b) *
            e^b - (1 + b) * x^b))/(e^b + x^b)^3
      },      
      inv = function(y, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            ((e^b * (-d + y))/(c - y))^(1/b) 
      },   
      expr.grad = expression(c + (d - c)/(1 + exp(b * (log(Cycles) - log(e))))), 
      inv.grad = expression(((e^b * (-d + Fluo))/(c - Fluo))^(1/b)), 
      parnames = c("b", "c", "d", "e"), 
      name = "l4",
      type = "four-parameter log-logistic"       
)

l3 <- list(
      expr = "Fluo ~ d/(1 + exp(b * (log(Cycles) - log(e))))",
      fct = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            d/(1 + exp(b * (log(x) - log(e))))
      },      
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/y2)
            lmFit <- rlm(logitTrans ~ log(x2))
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- exp(-coefVec[1]/b)
            ssVal <- as.numeric(c(b, d, e))
            names(ssVal) <- l3$parnames
            return(ssVal)
      },       
      d1 = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            -((b * d * e^b * x^(-1 + b))/(e^b + x^b)^2)
      },      
      d2 = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            -((b * d * e^b * x^(-2+b) * ((-1 + b) * e^b - (1 + b) * x^b))/(e^b + x^b)^3)
      },      
      inv = function(y, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            e * (-1 + (d/y))^(1/b)
      },   
      expr.grad = expression(d/(1 + exp(b * (log(Cycles) - log(e))))),
      inv.grad = expression(e * (-1 + (d/Fluo))^(1/b)),  
      parnames = c("b",  "d", "e"),
      name = "l3",
      type = "three-parameter log-logistic"
)

b7 <- list(
      expr = "Fluo ~ c + (k1 * Cycles) + (k2 * Cycles^2) + (d - c)/((1 + exp(b * (Cycles - e)))^f)",
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            c + (k1 * x) + (k2 * x^2) + (d - c)/((1 + exp(b * (x - e)))^f)
      },
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ x2)
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- -coefVec[1]/b
            f <- 1
            lmFit2 <- rlm(y2[1:10] ~ x2[1:10])
            k1 <- coef(lmFit2)[2]
            k2 <- - 0.01 * k1
            ssVal <- as.numeric(c(b, c, d, e, f, k1, k2))
            names(ssVal) <- b7$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            b * (c - d) * exp(b * (-e + x)) * (1 + exp(b * (-e + x)))^(-1 - f) * f + k1 + 2 * k2 * x
      },
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
            (1 + exp(b * (-e + x)))^(-2 - f) * (-b^2 * (c - d) * exp(b * (-e + x)) * f * (-1 + exp(b * (-e + x)) * f) +
            2 * (1 + exp(b * (-e + x)))^(2 + f) * k2)
      },
      inv = function(y, parm) {
            x <- 1:100
            fn <- function(x, parm) b7$fct(x, parm) - y 
            uniroot(fn, interval = c(1, 100), parm)$root
      },
      expr.grad = expression(c + (k1 * Cycles) + (k2 * Cycles^2) + (d - c)/((1 + exp(b * (Cycles - e)))^f)),
      inv.grad =  NULL,
      parnames = c("b", "c", "d", "e", "f", "k1", "k2"),
      name = "b7",
      type = "seven-parameter logistic"
)

b6 <- list(
      expr = "Fluo ~ c + (k * Cycles) + (d - c)/((1 + exp(b * (Cycles - e)))^f)", 
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            c + (k * x) + (d - c)/((1 + exp(b * (x - e)))^f)
      },
      ssFct = function (x, y) {
            d <- max(y) + 0.001                                     
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ x2)
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- -coefVec[1]/b
            f <- 1
            lmFit2 <- rlm(y2[1:10] ~ x2[1:10])
            k <- coef(lmFit2)[2]
            ssVal <- as.numeric(c(b, c, d, e, f, k))
            names(ssVal) <- b6$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            b * (c - d) * exp(b * (-e + x)) * (1 + exp(b * (-e + x)))^(-1-f) + k        
      },
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            -b^2 * (c - d) * exp(b * (-e + x)) * (1 + exp(b * (-e + x)))^(-2-f) * f * (-1 + exp(b * (-e + x)) * f)     
      },
      inv = function(y, parm) {
            x <- 1:100
            fn <- function(x, parm) b6$fct(x, parm) - y 
            uniroot(fn, interval = c(1, 100), parm)$root            
      },
      expr.grad = expression(c + (k * Cycles) + (d - c)/((1 + exp(b * (Cycles - e)))^f)),
      inv.grad =  NULL,
      parnames = c("b", "c", "d", "e", "f", "k"),
      name = "b6",
      type = "six-parameter logistic"
)

b5 <- list(
      expr = "Fluo ~ c + (d - c)/((1 + exp(b * (Cycles - e)))^f)", 
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            c + (d - c)/((1 + exp(b * (x - e)))^f)
      }, 
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ x2)
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- -coefVec[1]/b
            f <- 1
            ssVal <- as.numeric(c(b, c, d, e, f))
            names(ssVal) <- b5$parnames
            return(ssVal)
      },    
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            b * (c - d) * exp(b * (-e + x)) * (1 + exp(b * (-e + x)))^(-1 - f) * f   
      },   
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            -b^2 * (c - d) * exp(b * (-e + x)) * (1 + exp(b * (-e + x)))^(-2 - f) * 
            f * (-1 + exp(b * (-e + x)) * f)
      },    
      inv = function(y, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            e - ((log(1/(-1+((c-d)/(c-y))^(1/f))))/b) 
      },     
      expr.grad = expression(c + (d - c)/((1 + exp(b * (Cycles - e)))^f)),
      inv.grad = expression(e - ((log(1/(-1+((c-d)/(c-Fluo))^(1/f))))/b)), 
      parnames = c("b", "c", "d", "e", "f"),
      name = "b5",
      type = "five-parameter logistic"
)

b4 <- list(
      expr = "Fluo ~ c + (d - c)/(1 + exp(b * (Cycles - e)))", 
      fct = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            c + (d - c)/(1 + exp(b * (x - e)))
      },      
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            c <- min(y) - 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/(y2 - c))
            lmFit <- rlm(logitTrans ~ x2)
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- -coefVec[1]/b
            ssVal <- as.numeric(c(b, c, d, e))
            names(ssVal) <- b4$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            (b * (c - d) * exp(b * (e + x)))/(exp(b * e) + exp(b * x))^2
      }, 
      d2 = function(x, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            (b^2 * (c - d) * exp(b * (e + x)) * (exp(b * e) - exp(b * x)))/(exp(b * e) + exp(b * x))^3
      },  
      inv = function(y, parm) {
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            log((-d + y)/(c - y))/b + e
      }, 
      expr.grad = expression(c + (d - c)/(1 + exp(b * (Cycles - e)))), 
      inv.grad = expression(log((-d + Fluo)/(c - Fluo))/b + e),
      parnames = c("b", "c", "d", "e"),
      name = "b4",
      type = "four-parameter logistic"
)

b3 <- list(
      expr = "Fluo ~ d/(1 + exp(b * (Cycles - e)))",
      fct = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            d/(1 + exp(b * (x - e)))
      }, 
      ssFct = function (x, y) {
            d <- max(y) + 0.001
            x2 <- x[y > 0]
            y2 <- y[y > 0]
            logitTrans <- log((d - y2)/y2)
            lmFit <- rlm(logitTrans ~ x2)
            coefVec <- coef(lmFit)
            b <- coefVec[2]
            e <- -coefVec[1]/b
            ssVal <- as.numeric(c(b, d, e))
            names(ssVal) <- b3$parnames
            return(ssVal)
      },
      d1 = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            -(b * d * exp(b * (e + x)))/(exp(b * e) + exp(b * x))^2 
      },   
      d2 = function(x, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            (b^2 * d * exp(b * (e + x)) * (-exp(b * e) + exp(b * x)))/(exp(b * e) + exp(b * x))^3
      },              
      inv = function(y, parm) {
            b <- parm[1]
            d <- parm[2]
            e <- parm[3]
            log((d/y) - 1)/b + e
      },  
      expr.grad = expression(d/(1 + exp(b * (Cycles - e)))),
      inv.grad = expression(log((d/Fluo) - 1)/b + e),
      parnames = c("b",  "d", "e"),
      name = "b3",
      type = "three-parameter logistic"
)

expGrowth <- list(
  expr = "Fluo ~ a * exp(b * Cycles) + c", 
  fct = function(x, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    a * exp(b * x) + c
  },  
  ssFct = function (x, y) {  
    y[y == 0] <- 1E-12
    c <- 0.95 * min(y)     
    tempY <- log((y - c))    
    coefVec <- coef(rlm(tempY ~ x))
    a <- exp(coefVec[1])
    b <- coefVec[2] 
    ssVal <- c(a, b, c)
    names(ssVal) <- expGrowth$parnames             
    return(ssVal)
  },     
  d1 = function(x, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    a * b * exp(b * x)
  },     
  d2 = function(x, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    a * b^2 * exp(b * x)
  },        
  inv = function(y, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]               
    log((-c + y)/a)/b
  },       
  expr.grad = expression(a * exp(b * Cycles) + c), 
  inv.grad = expression(log((-c + Fluo)/a)/b),
  parnames = c("a", "b", "c"),       
  name = "expGrowth",
  type = "exponential growth model"
)

mak2 <- list(
  expr = "Fluo ~ mak2$fct(Cycles, c(D0, k, Fb))",
  fct = function(x, parm) {   
    D0 <- parm[1]
    k <- parm[2]
    if (k < 0.01) return(NA)
    Fb <- parm[3]      
    Fn <- vector(mode = "numeric", length = length(x))
    for (i in 1:length(x)) {
      if (i == 1) Fn[i] <- D0 + k * log(1 + (D0/k)) else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
    }
    Fn <- Fn + Fb    
    return(Fn)    
  },      
  ssFct = function(x, y) {
    sigDAT <- cbind(Cycles = x, Fluo = y)
    
    ## obtain cpD2 + offset and cut off all cycles beyond
    m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
    cpD2 <- efficiency(m, plot = FALSE)$cpD2 
    sigDAT <- sigDAT[1:floor(cpD2), ] 
      
    ## start estimates
    D0 <- 0.0001 
    k <- max(y, na.rm = TRUE)/10
    Fb <- min(y, na.rm = TRUE)    
    
    ssVal <- c(D0, k, Fb)
    names(ssVal) <- mak2$parnames    
    
    ## attach 'subset' attribute to result,
    ## so <pcrfit> knows which subset values to fit on.      
    attr(ssVal, "subset") <- 1:nrow(sigDAT)
    return(ssVal)    
  },
  d1 = function(x, parm) {           
  },
  d2 = function(x, parm) {            
  },
  inv = function(y, parm) {            
  },
  expr.grad = expression(mak2$fct(Cycles, c(D0, k, Fb))),
  inv.grad =  NULL,
  parnames = c("D0", "k", "Fb"),
  name = "mak2",
  type = "two-parameter mechanistic model"
)

mak2i <- list(
  expr = "Fluo ~ mak2i$fct(Cycles, c(D0, k, Fb))",
  fct = function(x, parm, y = NULL) {    
    D0 <- parm[1]
    k <- parm[2]
    if (k < 0.01) return(NA)
    Fb <- parm[3]             
    Fn <- vector(mode = "numeric", length = length(x))
    for (i in 1:length(x)) {
      if (i == 1) Fn[i] <- D0 + k * log(1 + (D0/k)) else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
    }
    Fn <- Fn + Fb
    if (is.null(y)) return(Fn) else return(y - Fn)    
  },      
  ssFct = function(x, y) {
    sigDAT <- cbind(Cycles = x, Fluo = y)
    
    ## obtain cpD2 + offset and cut off all cycles beyond
    m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
    cpD2 <- efficiency(m, plot = FALSE)$cpD2 
    sigDAT <- sigDAT[1:floor(cpD2), ]    
       
    ## make grid of start estimates
    D0.start <- 10^(-3:-12) 
    k.start <- seq(0.01, max(y, na.rm = TRUE)/10, length.out = 10)
    Fb.start <- min(y, na.rm = TRUE) 
        
    ## create grid of initial parameter values
    ## to optimize over
    START <- expand.grid(D0.start, k.start, Fb.start)
    
    ## initialize parameter matrix
    parMAT <- matrix(nrow = nrow(START), ncol = 4)
    colnames(parMAT) <- c("D0", "k", "Fb", "RSS")  
    
    ## nonlinear fitting of grid parameters
    for (i in 1:nrow(START)) {
      qpcR:::counter(i)
      PARM <- as.numeric(START[i, ]) 
      OUT <- try(nls.lm(par = PARM, fn = mak2i$fct, x = sigDAT[, 1], y = sigDAT[, 2], control = nls.lm.control(maxiter = 1000)), silent = TRUE)
      if (inherits(OUT, "try-error")) next
      RSS <- sum(OUT$fvec^2)
      parMAT[i, ] <- c(OUT$par, RSS)     
    }
    cat("\n")
    
    ## best value fit
    ORDER <- order(parMAT[, 4])
    parMAT <- parMAT[ORDER, ]
    parBEST <- parMAT[1, 1:3]
    names(parBEST) <- mak2i$parnames
    
    ## attach 'subset' attribute to result,
    ## so <pcrfit> knows which subset values to fit on.      
    attr(parBEST, "subset") <- 1:nrow(sigDAT)
    return(parBEST)
  },
  d1 = function(x, parm) {           
  },
  d2 = function(x, parm) {            
  },
  inv = function(y, parm) {            
  },
  expr.grad = expression(mak2i$fct(Cycles, c(D0, k, Fb))),
  inv.grad =  NULL,
  parnames = c("D0", "k", "Fb"),
  name = "mak2i",
  type = "(grid-searched) two-parameter mechanistic model"
)

mak3 <- list(
  expr = "Fluo ~ mak3$fct(Cycles, c(D0, k, Fb, slope))",
  fct = function(x, parm) {    
    D0 <- parm[1]
    k <- parm[2]
    if (k < 0.01) return(NA)
    Fb <- parm[3] 
    slope <- parm[4]
    Fn <- vector(mode = "numeric", length = length(x))
    for (i in 1:length(x)) {
      if (i == 1) Fn[i] <- D0 + k * log(1 + (D0/k)) else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
    }
    Fn <- Fn + Fb
    Fn <- Fn + (slope * (1:length(Fn)))
    
    return(Fn)    
  },      
  ssFct = function(x, y) {
    sigDAT <- cbind(Cycles = x, Fluo = y)
    
    ## obtain cpD2 + offset and cut off all cycles beyond
    m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
    cpD2 <- efficiency(m, plot = FALSE)$cpD2 
    sigDAT <- sigDAT[1:floor(cpD2), ]    
    
    ## start estimates
    D0 <- 0.0001 
    k <- max(y, na.rm = TRUE)/10
    Fb <- min(y, na.rm = TRUE)
    slope <- coef(lm(sigDAT[1:8, 2] ~ I(1:8)))[2]
    
    ssVal <- c(D0, k, Fb, slope)
    names(ssVal) <- mak3$parnames
    
    ## attach 'subset' attribute to result,
    ## so <pcrfit> knows which subset values to fit on.      
    attr(ssVal, "subset") <- 1:nrow(sigDAT)
    return(ssVal)
  },
  d1 = function(x, parm) {           
  },
  d2 = function(x, parm) {            
  },
  inv = function(y, parm) {            
  },
  expr.grad = expression(mak3$fct(Cycles, c(D0, k, Fb, slope))),
  inv.grad =  NULL,
  parnames = c("D0", "k", "Fb", "slope"),
  name = "mak3",
  type = "three-parameter mechanistic model"
)

mak3i <- list(
  expr = "Fluo ~ mak3i$fct(Cycles, c(D0, k, Fb, slope))",
  fct = function(x, parm, y = NULL) {    
    D0 <- parm[1]
    k <- parm[2]
    if (k < 0.01) return(NA)
    Fb <- parm[3] 
    slope <- parm[4]
    Fn <- vector(mode = "numeric", length = length(x))
    for (i in 1:length(x)) {
      if (i == 1) Fn[i] <- D0 + k * log(1 + (D0/k)) else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
    }
    Fn <- Fn + Fb
    Fn <- Fn + (slope * (1:length(Fn)))
    
    if (is.null(y)) return(Fn) else return(y - Fn)    
  },      
  ssFct = function(x, y) {
    sigDAT <- cbind(Cycles = x, Fluo = y)
    
    ## obtain cpD2 + offset and cut off all cycles beyond
    m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
    cpD2 <- efficiency(m, plot = FALSE)$cpD2 
    sigDAT <- sigDAT[1:floor(cpD2), ]    
    
    ## make grid of start estimates
    D0.start <- 10^(-3:-12) 
    k.start <- seq(0.01, max(y, na.rm = TRUE)/10, length.out = 10)
    Fb.start <- min(y, na.rm = TRUE) 
    slope <- coef(lm(sigDAT[1:8, 2] ~ I(1:8)))[2]
    
    ## create grid of initial parameter values
    ## to optimize over
    START <- expand.grid(D0.start, k.start, Fb.start, slope)
    
    ## initialize parameter matrix
    parMAT <- matrix(nrow = nrow(START), ncol = 5)
    colnames(parMAT) <- c("D0", "k", "Fb", "slope", "RSS")  
    
    ## nonlinear fitting of grid parameters
    for (i in 1:nrow(START)) {
      qpcR:::counter(i)
      PARM <- as.numeric(START[i, ]) 
      OUT <- try(nls.lm(par = PARM, fn = mak3i$fct, x = sigDAT[, 1], y = sigDAT[, 2], control = nls.lm.control(maxiter = 1000)), silent = TRUE)
      if (inherits(OUT, "try-error")) next
      RSS <- sum(OUT$fvec^2)
      parMAT[i, ] <- c(OUT$par, RSS)     
    }
    cat("\n")
    
    ## best value fit
    ORDER <- order(parMAT[, 5])
    parMAT <- parMAT[ORDER, ]
    parBEST <- parMAT[1, 1:4]
    names(parBEST) <- mak3i$parnames
    
    ## attach 'subset' attribute to result,
    ## so <pcrfit> knows which subset values to fit on.      
    attr(parBEST, "subset") <- 1:nrow(sigDAT)
    return(parBEST)
  },
  d1 = function(x, parm) {           
  },
  d2 = function(x, parm) {            
  },
  inv = function(y, parm) {            
  },
  expr.grad = expression(mak3i$fct(Cycles, c(D0, k, Fb, slope))),
  inv.grad =  NULL,
  parnames = c("D0", "k", "Fb", "slope"),
  name = "mak3i",
  type = "(grid-searched) three-parameter mechanistic model"
)

lin2 <- list(
  expr = "Fluo ~ eta * log(exp(a1 * (Cycles - tau)/eta) + exp(a2 * (Cycles - tau)/eta)) + c", 
  fct = function(x, parm) {
    c <- parm[1]
    eta <- parm[2]
    tau <- parm[3]
    a1 <- parm[4]
    a2 <- parm[5]
    eta * log(exp(a1 * (x - tau)/eta) + exp(a2 * (x - tau)/eta)) + c
  },  
  ssFct = function (x, y) {
    sigDAT <- cbind(Cycles = x, Fluo = y)
    m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)
    cpD2 <- efficiency(m, plot = FALSE)$cpD2
    cpD2 <- floor(cpD2)
    c <- min(y) - 0.001
    sub1 <- 1:(cpD2 - 5)
    fit1 <- lm(y[sub1] ~ x[sub1])  
    a1 <- coef(fit1)[2]
    sub2 <- cpD2:(cpD2 + 3)
    fit2 <- lm(y[sub2] ~ x[sub2])    
    a2 <- coef(fit2)[2]
    tau <- cpD2
    eta <- a2 - a1
    ssVal <- as.numeric(c(c, eta, tau, a1, a2))
    names(ssVal) <- lin2$parnames
    attr(ssVal, "subset") <- 1:(cpD2 + 3)
    return(ssVal)
  },     
  d1 = function(x, parm) {
    c <- parm[1]
    eta <- parm[2]
    tau <- parm[3]
    a1 <- parm[4]
    a2 <- parm[5]
    eta * ((exp(a1 * (x - tau)/eta) * (a1/eta) + 
      exp(a2 * (x - tau)/eta) * (a2/eta))/(exp(a1 * (x - tau)/eta) + 
      exp(a2 * (x - tau)/eta)))    
  },        
  d2 = function(x, parm) {
    c <- parm[1]
    eta <- parm[2]
    tau <- parm[3]
    a1 <- parm[4]
    a2 <- parm[5]
    eta * ((exp(a1 * (x - tau)/eta) * (a1/eta) * (a1/eta) + 
      exp(a2 * (x - tau)/eta) * (a2/eta) * (a2/eta))/(exp(a1 * 
      (x - tau)/eta) + exp(a2 * (x - tau)/eta)) - (exp(a1 * 
      (x - tau)/eta) * (a1/eta) + exp(a2 * (x - tau)/eta) * 
      (a2/eta)) * (exp(a1 * (x - tau)/eta) * (a1/eta) + exp(a2 * 
      (x - tau)/eta) * (a2/eta))/(exp(a1 * (x - tau)/eta) + 
      exp(a2 * (x - tau)/eta))^2)    
  },  
  inv = function(y, parm) {
    x <- 1:100
    fn <- function(x, parm) lin2$fct(x, parm) - y
    uniroot(fn, interval = c(1, 100), parm)$root
  },
  expr.grad = NULL,
  inv.grad =  NULL,      
  parnames = c("c", "eta", "tau", "a1", "a2"),
  name = "lin2",
  type = "bilinear model"
)

cm3 <- list(
  expr = "Fluo ~ cm3$fct(Cycles, c(D0, max, Kd, Fb))",
  fct = function(x, parm) {    
    D0 <- parm[1]
    if (D0 <= 0) return(NA)
    max <- parm[2]
    Kd <- parm[3] 
    if (Kd <= 0) return(NA)
    Fb <- parm[4]
    Fn <- vector(mode = "numeric", length = length(x))
    for (i in 1:length(x)) {
      if (i == 1) Fn[i] <- D0 * (1 + ((max - D0)/max) - (D0/(Kd + D0))) else Fn[i] <- Fn[i - 1] * (1 + ((max - Fn[i - 1])/max) - (Fn[i - 1]/(Kd + Fn[i - 1])))  
    }
    Fn <- Fn + Fb
    return(Fn)    
  },      
  ssFct = function(x, y) {
    ## start estimates
    D0 <- 0.00001 
    max <- 3 * max(y, na.rm = TRUE)
    Kd <- 0.5 * max(y, na.rm = TRUE)
    LM <- lm(y[1:8] ~ x[1:8]) 
    Fb <- coef(LM)[2] + 0.01
    ssVal <- c(D0, max, Kd, Fb)
    names(ssVal) <- cm3$parnames    
    return(ssVal)
  },
  d1 = function(x, parm) {           
  },
  d2 = function(x, parm) {            
  },
  inv = function(y, parm) {
    x <- 1:100
    fn <- function(x, parm) cm3$fct(x, parm) - y
    uniroot(fn, interval = c(1, 100), parm)$root
  },
  expr.grad = expression(cm3$fct(Cycles, c(D0, max, Kd, Fb))),
  inv.grad =  NULL,
  parnames = c("D0", "max", "Kd", "Fb"),
  name = "cm3",
  type = "three-parameter mechanistic model"
)

linexp <- list(
  expr = "Fluo ~ a * exp(b * Cycles) + (k * Cycles) + c",
  fct = function (x, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    k <- parm[4]
    a * exp(b * x) + (k * x) + c
  },
  ssFct = function (x, y) {
    DATA <- data.frame(Cycles = x, Fluo = y)
    FIT <- pcrfit(DATA, 1, 2, expGrowth, verbose = FALSE)
    COEF <- coef(FIT)
    lmFit <- rlm(y[1:10] ~ x[1:10])
    k <- coef(lmFit)[2]
    ssVal <- c(COEF, k)
    names(ssVal) <- linexp$parnames
    return(ssVal)
  },
  d1 = function (x, parm) {
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    k <- parm[4]
    a * (exp(b * x) * b) + k
  },
  d2 = function (x, parm) { 
    a <- parm[1]
    b <- parm[2]
    c <- parm[3]
    k <- parm[4]
    a * (exp(b * x) * b * b)
  },
  inv = function (y, parm) {
    x <- 1:100
    fn <- function(x, parm) linexp$fct(x, parm) - y
    uniroot(fn, interval = c(1, 100), parm)$root
  },
  expr.grad = expression(a * exp(b * Cycles) + (k * Cycles) + c),
  inv.grad = NULL,
  parnames = c("a", "b", "c", "k"),
  name = "linexp",
  type = "linear-exponential growth model"
)

spl3 <- list(
  expr = "Fluo ~ spl3$fct(Cycles, c(spar), Fluo)",
  fct = function(x, parm, yvec = NULL) { 
    spar <- parm[1]   
    if (length(x) != length(yvec)) {     
      SPL <- smooth.spline(1:length(yvec), yvec, spar = spar) 
      PRED <- predict(SPL, x)
      return(PRED$y)
    }
    SPL <- smooth.spline(x, yvec, spar = spar)
    Fn <- SPL$y   
    return(Fn)   
  },     
  ssFct = function(x, y) {
    ## start estimates
    spar <- 1
    assign("yvec", y, envir = .GlobalEnv)
    ssVal <- c(spar) 
    names(ssVal) <- spl3$parnames   
    return(ssVal)
  },  
  d1 = function(x, parm) {          
  },
  d2 = function(x, parm) {           
  },
  inv = function(y, parm, yvec) {
    x <- 1:100
    fn <- function(x, parm, yvec) spl3$fct(x, parm, yvec) - y    
    uniroot(fn, interval = c(1, 100), parm, yvec)$root
  },
  expr.grad = expression(spl3$fct(Cycles, c(spar))),
  inv.grad =  NULL,
  parnames = c("spar"),
  name = "spl3",
  type = "cubic smoothing spline"
)

########################
save.image(file = "c:\\temp\\qpcR_functions.rda")


