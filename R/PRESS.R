PRESS <- function(model)
{
      VARS <- all.vars(model$call$formula)
      CALL <- as.list(model$call)

      if (!is.null(as.list(model$call)$data)) {
            if (!is.data.frame(as.list(model$call)$data)) {
                  MATCH <- which(as.list(model$call)$data == all.vars(model$call))
                  DATA <- get(noquote(all.vars(model$call)[MATCH]))
                  suppressMessages(attach(DATA))
            } else {
                  DATA <- as.list(model$call)$data
                  suppressMessages(attach(DATA))
            }
      }
      
      DATA <- as.data.frame(sapply(VARS, function(a) get(noquote(a))))
      PRESS.res <- NULL

      for (i in 1:nrow(DATA)) {
            NEWCALL <- CALL
            NEWCALL$data <- DATA[-i, ]
            NEWMOD <- eval(as.call(NEWCALL))
            NEWPRED <- as.data.frame(DATA[i, -1])
            colnames(NEWPRED) <- VARS[-1]
            y.hat <- predict(NEWMOD, NEWPRED)[1]
            PRESS.res[i] <- DATA[i, 1] - y.hat
      }
      
      return(list(stat = sum(PRESS.res^2), residuals = PRESS.res))
}