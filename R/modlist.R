modlist <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
opt = FALSE, 
norm = FALSE,
backsub = NULL, 
opt.method =  "LM",
nls.method = "port",
sig.level = 0.05, 
crit = "ftest",
...
)
{
  modList <- NULL
  counter <- 1   
  
  if (class(x)[1] == "pcrfit") {
    model <- x$MODEL
    x <- x$DATA            
  }                         
  
  if (is.null(fluo)) fluo <- 2:ncol(x)        
    
  if (!is.null(backsub) && !is.numeric(backsub)) 
      stop("'backsub' must be either NULL or a numeric sequence!")     

  for (i in fluo) {
    Cycles <- x[, cyc]
    Fluo  <- x[, i]
    compl <- complete.cases(Fluo)
    Cycles <- Cycles[compl]
    Fluo <- Fluo[compl]
  
    if (norm) {
      Fluo <- Fluo - min(Fluo, na.rm = TRUE)
      Fluo <- Fluo/max(Fluo, na.rm = TRUE)
    }
        
    if (!is.null(backsub)) {
      back <- mean(Fluo[backsub], na.rm = TRUE)
      Fluo <- Fluo - back
    }                    
        
    DATA <- data.frame(Cycles = Cycles, Fluo = Fluo)     
    flush.console()
    cat("Making model for ", colnames(x[i]), " (", model$name, ")", sep= "")
                
    fitObj <- try(pcrfit(DATA, 1, 2, model, opt.method = opt.method, nls.method = nls.method, ...), silent = TRUE)
    
    if (inherits(fitObj, "try-error")) {
      cat(" => ", colnames(x[i]), " gave a fitting error!", sep = "")
      fitObj$DATA <- DATA
      class(fitObj) <- "pcrfit"       
    }
      
    if (opt) {
 	    fitObj2 <- try(mselect(fitObj, verbose = FALSE, sig.level = sig.level, crit = crit, ...), silent = TRUE)             
      if (inherits(fitObj2, "try-error")) {
        fitObj <- fitObj
        cat(" => ", colnames(x[i]), " gave a model selection error!", sep = "")
      } else {
        fitObj <- fitObj2
        cat(" => ", fitObj$MODEL$name, sep = "")
      }
    }   
    
    cat("\n")
           
    modList[[counter]] <- fitObj
    modList[[counter]]$names <- colnames(x[i])   
    counter <- counter + 1       
  }
    
  class(modList) <- c("modlist", "pcrfit")
  invisible(modList)
}

