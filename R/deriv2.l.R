deriv2.l <-
function (x, parm, ffun) 
    {
     if (ffun == "l5"){fixed <- c(NA, NA, NA, NA, NA)}
     if (ffun == "l4"){fixed <- c(NA, NA, NA, NA, 1)}
     if (ffun == "l3"){fixed <- c(NA, 0, NA, NA, 1)}
     notFixed <- is.na(fixed)
     parmVec <- rep(0, numParm <- 5 )
     parmVec[!notFixed] <- fixed[!notFixed]
     parmMat <- matrix(parmVec, nrow(parm), numParm <- 5 , byrow = TRUE)
     parmMat[, notFixed] <- parm
    .expr1 <-  parmMat[,3] -  parmMat[,2]
    .expr3 <- x/parmMat[,4]
    .expr5 <- 1 + .expr3^parmMat[,1]
    .expr6 <- parmMat[,5] - 1
    .expr9 <-  parmMat[,1]/ parmMat[,4]
    .expr10 <-  parmMat[,5] * .expr5^.expr6 * .expr9
    .expr11 <-  parmMat[,1] - 1
    .expr12 <- .expr3^.expr11
    .expr14 <- -.expr1 * (.expr10 * .expr12)
    .expr15 <- 2 *  parmMat[,5]
    .expr16 <- .expr5^.expr15
    .expr20 <- 1/ parmMat[,4]
    .expr22 <- .expr12 * ( parmMat[,1] * .expr20)
    .value <- .expr14/.expr16
    -(.expr1 * (parmMat[,5] * (.expr5^(.expr6 - 1) * (.expr6 * 
        .expr22)) * .expr9 * .expr12 + .expr10 * (.expr3^(.expr11 - 
        1) * (.expr11 * .expr20)))/.expr16 + .expr14 * (.expr5^(.expr15 - 
        1) * (.expr15 * .expr22))/.expr16^2)
     }
