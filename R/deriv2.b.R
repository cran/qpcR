deriv2.b <-
function (x, parm, ffun) 
    {
     if (ffun == "b5") {fixed <- c(NA, NA, NA, NA, NA)}
     if (ffun == "b4") {fixed <- c(NA, NA, NA, NA, 1)}
     if (ffun == "b3") {fixed <- c(NA, 0, NA, NA, 1)}
     notFixed <- is.na(fixed)
     parmVec <- rep(0, numParm <- 5 )
     parmVec[!notFixed] <- fixed[!notFixed]
     parmMat <- matrix(parmVec, nrow(parm), numParm <- 5 , byrow = TRUE)
     parmMat[, notFixed] <- parm
    .expr1 <- parmMat[,3] - parmMat[,2]
    .expr5 <- exp(parmMat[,1] * (x - parmMat[,4]))
    .expr7 <- -.expr1 * .expr5 * parmMat[,1]
    .expr8 <- 1 + .expr5
    .expr9 <- .expr8^2
    .expr11 <- .expr5 * parmMat[,1]
    .value <- .expr7/.expr9
    -(.expr1 * .expr11 * parmMat[,1]/.expr9 + .expr7 * (2 * (.expr11 * .expr8))/.expr9^2)
     }

