modlist <- function(
x, 
cyc = 1, 
fluo = NULL, 
model = l4, 
opt = FALSE, 
norm = FALSE,
backsub = NULL, 
opt.method =  rep("Nelder", 5),
nls.method = "port",
sig.level = 0.05, 
crit = "ftest",
...
)
{
  modList <- NULL
  counter <- 1                   
  
  if (is.null(fluo)) fluo <- 2:ncol(x)        
    
  if (!is.null(backsub) && !is.numeric(backsub)) 
      stop("'backsub' must be either NULL or a numeric sequence!") 
    
  for (i in fluo) {
    Cycles <- x[, cyc]
    Fluo  <- x[, i]
               
    if (norm) {
      Fluo <- Fluo - min(Fluo, na.rm = TRUE)
      Fluo <- Fluo/max(Fluo, na.rm = TRUE)
    }
        
    if (!is.null(backsub)) {
      back <- mean(Fluo[backsub], na.rm = TRUE)
      Fluo <- Fluo - back
    }                    
        
    DATA <- cbind(Cycles, Fluo)  
    flush.console()
    cat("Making model for ", colnames(x[i]), " (", model$name, ")\n", sep= "") 
                
    fitObj <- pcrfit(DATA, 1, 2, model, opt.method = opt.method, nls.method = nls.method, ...)          
        
    if (opt) {
 	    fitObj <- try(mselect(fitObj, verbose = FALSE, sig.level = sig.level, crit = crit, ...), silent = TRUE)             
    }                     
        
    modList[[counter]] <- fitObj
    modList[[counter]]$names <- colnames(x[i])   
    counter <- counter + 1       
  }
    
  class(modList) <- c("modlist", "pcrfit")
  invisible(modList)
}
