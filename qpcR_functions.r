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
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
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
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]             
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
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k1 <- parm[6]
            k2 <- parm[7]
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
            b <- parm[1]
            c <- parm[2]
            d <- parm[3]
            e <- parm[4]
            f <- parm[5]
            k <- parm[6]
            
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
            plateau <- 0.95 * min(y)
            span <- max(y) - plateau
            tempY <- log((y - plateau))
            coefVec <- coef(rlm(tempY ~ x))
            span <- exp(coefVec[1])
            K <- coefVec[2] 
            ssVal <- c(span, K, plateau)
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

mak3 <- list(
      expr = "Fluo ~ mak3$fct(Cycles, ssVal)",
      fct = function(x, parm) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]
            slope <- parm[4]
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + (slope * (1:length(Fn)) + Fb)
            return(Fn)
      },         
      fct_ssFct = function(parm, x, y, method) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]
            slope <- parm[4]
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + (slope * (1:length(Fn)) + Fb)
            if (method == "LM") res <- y - Fn else res <- sum(y - Fn)^2
            return(res)            
      },
      ssFct = function(x, y) {
            ## get parMAKs parameters from global environment or use
            ## standard ones
            if (exists("parMAKs", envir = .GlobalEnv)) {
              PARS <- get("parMAKs", envir = .GlobalEnv)
            } else PARS <- parMAK()            
            
            D2.offset <- PARS$SS.offset 
            method <- PARS$SS.method 
            cutter <- PARS$SS.deriv          
            cutter <- match.arg(cutter, c("sigfit", "spline"))
            
            # select cut-off method
            sigDAT <- cbind(Cycles = x, Fluo = y)
            if (cutter == "sigfit") {
              t1 <- supsmu(x = 1:length(y), y = y, span = 0.1)$y
              t2 <- diff(t1)
              t3 <- diff(t2)
              cpD2 <- which.max(t3) + D2.offset
            } else {               
              m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
              cpD2 <- efficiency(m, plot = FALSE)$cpD2 + D2.offset 
            }              
            # cut off all cycles beyond...
            sigDAT <- sigDAT[1:floor(cpD2), ]
            # exchange 0 with small value
            sigDAT[,2][sigDAT[, 2] == 0] <- 1E-6
            # exponential fit for Fb and D0 start estimates
            m2 <- pcrfit(sigDAT, 1, 2, expGrowth, verbose = FALSE)
            # make grid of start estimates
            D0.start <- coef(m2)[1] * 10^(-4:1)
            k.start <- seq(0.1, 3, by = 0.3)
            Fb.start <- coef(m2)[3]
            slope.start <- coef(lm(sigDAT[1:5, 2] ~ I(1:5)))[2]
            ### create grid of initial parameter values
            ### to optimize over
            START <- expand.grid(D0.start, k.start, Fb.start, slope.start)
            ### initialize parameter matrix
            parMAT <- matrix(nrow = nrow(START), ncol = 5)
            colnames(parMAT) <- c("D0", "k", "Fb", "slope", "RSS")              
            ### nonlinear fitting
            for (i in 1:nrow(START)) {
              PARM <- as.numeric(START[i, ])
              if (method == "LM") OUT <- try(nls.lm(PARM, mak3$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = nls.lm.control(maxiter = 1000)), silent = TRUE)
              else OUT <- try(optim(PARM, mak3$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = list(maxit = 1000)), silent = TRUE)
              if (inherits(OUT, "try-error")) next
              RSS <- if(method == "LM") sum(OUT$fvec^2) else OUT$value
              parMAT[i, ] <- c(OUT$par, RSS)
              qpcR:::counter(i)
            }
            cat("\n")
            ### best value fit
            ORDER <- order(parMAT[, 5])
            parMAT <- parMAT[ORDER, ]
            parBEST <- parMAT[1, 1:4]
            
            names(parBEST) <- mak3$parnames
                         
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
      expr.grad = expression(mak3$fct(Cycles, ssVal)),
      inv.grad =  NULL,
      parnames = c("D0", "k", "Fb", "slope"),
      name = "mak3",
      type = "three-parameter mechanistic model"
)

mak2 <- list(
      expr = "Fluo ~ mak2$fct(Cycles, ssVal)",
      fct = function(x, parm) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]             
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + Fb
            return(Fn)
      },         
      fct_ssFct = function(parm, x, y, method) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]              
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + Fb
            if (method == "LM") res <- y - Fn else res <- sum(y - Fn)^2
            return(res)            
      },
      ssFct = function(x, y) {
            ## get parMAKs parameters from global environment or use
            ## standard ones
            if (exists("parMAKs", envir = .GlobalEnv)) {
              PARS <- get("parMAKs", envir = .GlobalEnv)
            } else PARS <- parMAK()
                       
            D2.offset <- PARS$SS.offset 
            method <- PARS$SS.method 
            cutter <- PARS$SS.deriv          
            cutter <- match.arg(cutter, c("sigfit", "spline"))
                        
            # select cut-off method
            sigDAT <- cbind(Cycles = x, Fluo = y)
                        
            if (cutter == "sigfit") {
              t1 <- supsmu(x = 1:length(y), y = y, span = 0.1)$y
              t2 <- diff(t1)
              t3 <- diff(t2)
              cpD2 <- which.max(t3) + D2.offset
            } else {             
              m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
              cpD2 <- efficiency(m, plot = FALSE)$cpD2 + D2.offset 
            }
                     
            # cut off all cycles beyond...
            sigDAT <- sigDAT[1:floor(cpD2), ]
            # exchange 0 with small value
            sigDAT[,2][sigDAT[, 2] == 0] <- 1E-6
            # exponential fit for Fb and D0 start estimates
            m2 <- pcrfit(sigDAT, 1, 2, expGrowth, verbose = FALSE)
            # make grid of start estimates
            D0.start <- coef(m2)[1] * 10^(-4:1)
            k.start <- seq(0.1, 3, by = 0.3)
            Fb.start <- coef(m2)[3]             
            ### create grid of initial parameter values
            ### to optimize over
            START <- expand.grid(D0.start, k.start, Fb.start)
            ### initialize parameter matrix
            parMAT <- matrix(nrow = nrow(START), ncol = 4)
            colnames(parMAT) <- c("D0", "k", "Fb", "RSS")              
            ### nonlinear fitting
            for (i in 1:nrow(START)) {
              PARM <- as.numeric(START[i, ])            
              if (method == "LM") OUT <- try(nls.lm(PARM, mak2$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = nls.lm.control(maxiter = 1000)), silent = TRUE)
              else OUT <- try(optim(PARM, mak2$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = list(maxit = 1000)), silent = TRUE)
              if (inherits(OUT, "try-error")) next
              RSS <- if(method == "LM") sum(OUT$fvec^2) else OUT$value
              parMAT[i, ] <- c(OUT$par, RSS)
              qpcR:::counter(i)
            }
            cat("\n")
            ### best value fit
            ORDER <- order(parMAT[, 4])
            parMAT <- parMAT[ORDER, ]
            parBEST <- parMAT[1, 1:3]
            names(parBEST) <- mak2$parnames
                         
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
      expr.grad = expression(mak2$fct(Cycles, ssVal)),
      inv.grad =  NULL,
      parnames = c("D0", "k", "Fb"),
      name = "mak2",
      type = "two-parameter mechanistic model"
)

chag <- list(
      expr = "Fluo ~ chag$fct(Cycles, ssVal)",
      fct = function(x, parm) {
            d0 <- parm[1]
            a <- parm[2]
            b <- parm[3]            
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- 1 - (1 - Fn[i-1]) * ((1 - b * Fn[i-1])/(1 + (a - 2) * b * Fn[i-1]))^(1/(a-1))
            }
            return(Fn)
      },         
      fct_ssFct = function(parm, x, y, method) {
            d0 <- parm[1]
            a <- parm[2]
            b <- parm[3]
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                 if (i == 1) Fn[i] <- d0 else Fn[i] <- 1 - (1 - Fn[i-1]) * ((1 - b * Fn[i-1])/(1 + (a - 2) * b * Fn[i-1]))^(1/(a-1))
            }            
            if (method == "LM") res <- y - Fn else res <- sum(y - Fn)^2
            return(res)            
      },
      ssFct = function(x, y) {   
            ## get parMAKs parameters from global environment or use
            ## standard ones
            if (exists("parMAKs", envir = .GlobalEnv)) {
              PARS <- get("parMAKs", envir = .GlobalEnv)
            } else PARS <- parMAK()             
           
            method <- PARS$SS.method              
        
            ### baseline subtraction of cycles 3-8
            LM <- lm(y[3:8] ~ x[3:8])
            y <- y - coef(LM)[1]        
            ### normalize fluorescence values within [0, 1]
            y <- qpcR:::rescale(y, 0, 1)           
            # make grid of start estimates
            D0.start <- 10^(-3:-12)
            a.start <- seq(1, 10, by = 3)
            b.start <- seq(0.5, 2, by = 0.5)            
            ### create grid of initial parameter values
            ### to optimize over
            START <- expand.grid(D0.start, a.start, b.start)            
            ### initialize parameter matrix
            parMAT <- matrix(nrow = nrow(START), ncol = 4)
            colnames(parMAT) <- c("D0", "a", "b", "RSS")              
            ### nonlinear fitting
            for (i in 1:nrow(START)) {
              PARM <- as.numeric(START[i, ])
              if (method == "LM") OUT <- try(nls.lm(PARM, chag$fct_ssFct, x = x, y = y, method = method, control = nls.lm.control(maxiter = 1000)), silent = TRUE)
              else OUT <- try(optim(PARM, chag$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = list(maxit = 1000)), silent = TRUE)
              if (inherits(OUT, "try-error")) next
              RSS <- sum(OUT$fvec^2)
              parMAT[i, ] <- c(OUT$par, RSS)
              qpcR:::counter(i)
            }            
            cat("\n")
                        
            ### best value fit
            ORDER <- order(parMAT[, 4])
            parMAT <- parMAT[ORDER, ]
            parBEST <- parMAT[1, 1:3]            
            
            names(parBEST) <- chag$parnames       
            
            ## attach 'subset' attribute to result,
            ## so <pcrfit> knows which subset values to fit on.      
            attr(parBEST, "scale") <- c(0, 1)
            return(parBEST)
      },
      d1 = function(x, parm) {                 
      },
      d2 = function(x, parm) {            
      },
      inv = function(y, parm) {            
      },
      expr.grad = expression(chag$fct(Cycles, ssVal)),
      inv.grad =  NULL,
      parnames = c("D0", "a", "b"),
      name = "chag",
      type = "Chagovetz's iterative map"
)

mak3n <- list(
      expr = "Fluo ~ mak3n$fct(Cycles, ssVal)",
      fct = function(x, parm) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]
            slope <- parm[4]
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + (slope * (1:length(Fn)) + Fb)
            return(Fn)
      },         
      fct_ssFct = function(parm, x, y, method = method) {
            d0 <- parm[1]
            k <- parm[2]
            Fb <- parm[3]
            slope <- parm[4]
            Fn <- vector(mode = "numeric", length = length(x))
            for (i in 1:length(x)) {
                if (i == 1) Fn[i] <- d0 else Fn[i] <- Fn[i-1] + k * log(1 + (Fn[i-1]/k))
            }
            Fn <- Fn + (slope * (1:length(Fn)) + Fb)
            if (method == "LM") res <- y - Fn else res <- sum(y - Fn)^2
            return(res)            
      },
      ssFct = function(x, y) {
            
            ## get parMAKs parameters from global environment or use
            ## standard ones
            if (exists("parMAKs", envir = .GlobalEnv)) {
              PARS <- get("parMAKs", envir = .GlobalEnv)
            } else PARS <- parMAK()
                       
            D2.offset <- PARS$SS.offset 
            method <- PARS$SS.method 
            cutter <- PARS$SS.deriv          
            cutter <- match.arg(cutter, c("sigfit", "spline"))
                
            # select cut-off method
            y <- qpcR:::rescale(y, 0, 1)
            sigDAT <- cbind(Cycles = x, Fluo = y)
            if (cutter == "sigfit") {
              t1 <- supsmu(x = 1:length(y), y = y, span = 0.1)$y
              t2 <- diff(t1)
              t3 <- diff(t2)
              cpD2 <- which.max(t3) + D2.offset
            } else {               
              m <- pcrfit(sigDAT, 1, 2, l4, verbose = FALSE)             
              cpD2 <- efficiency(m, plot = FALSE)$cpD2 + D2.offset 
            }              
            # cut off all cycles beyond...
            sigDAT <- sigDAT[1:floor(cpD2), ]
            # exchange 0 with small value
            sigDAT[,2][sigDAT[, 2] == 0] <- 1E-6
            # exponential fit for Fb and D0 start estimates
            #m2 <- pcrfit(sigDAT, 1, 2, expGrowth, verbose = FALSE)
            # make grid of start estimates
            D0.start <- 10^-(4:12)
            k.start <- 0.01
            Fb.start <- 0.001
            slope.start <- 0.001
            ### create grid of initial parameter values
            ### to optimize over
            START <- expand.grid(D0.start, k.start, Fb.start, slope.start)            
            ### initialize parameter matrix
            parMAT <- matrix(nrow = nrow(START), ncol = 5)
            colnames(parMAT) <- c("D0", "k", "Fb", "slope", "RSS")              
            ### nonlinear fitting
            for (i in 1:nrow(START)) {
              PARM <- as.numeric(START[i, ])
              if (method == "LM") OUT <- try(nls.lm(PARM, mak3n$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = nls.lm.control(maxiter = 1000)), silent = TRUE)
              else OUT <- try(optim(PARM, mak3n$fct_ssFct, x = sigDAT[, 1], y = sigDAT[, 2], method = method, control = list(maxit = 1000)), silent = TRUE)
              if (inherits(OUT, "try-error")) next
              RSS <- if(method == "LM") sum(OUT$fvec^2) else OUT$value
              parMAT[i, ] <- c(OUT$par, RSS)
              qpcR:::counter(i)
            }
            cat("\n")
            ### best value fit
            ORDER <- order(parMAT[, 5])
            parMAT <- parMAT[ORDER, ]
            parBEST <- parMAT[1, 1:4]
            
            names(parBEST) <- mak3n$parnames
                         
            ## attach 'subset' attribute to result,
            ## so <pcrfit> knows which subset values to fit on.
            attr(parBEST, "scale") <- c(0, 1)
            attr(parBEST, "subset") <- 1:nrow(sigDAT)
            return(parBEST)
      },
      d1 = function(x, parm) {                 
      },
      d2 = function(x, parm) {            
      },
      inv = function(y, parm) {            
      },
      expr.grad = expression(mak3n$fct(Cycles, ssVal)),
      inv.grad =  NULL,
      parnames = c("D0", "k", "Fb", "slope"),
      name = "mak3",
      type = "three-parameter mechanistic model"
)

########################
save.image(file = "c:\\temp\\qpcR_functions.rda")


