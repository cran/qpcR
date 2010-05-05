counter <- function(i) {
  if (i %% 10 == 0) cat(i) else cat(".")
  if (i %% 50 == 0) cat("\n")
  flush.console()
}

