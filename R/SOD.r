SOD <- function(
object,    
remove = FALSE,
verbose = TRUE, 
...
)
{
 CLASS <- class(object)

 if (!(CLASS[1] %in% c("pcrfit", "modlist", "replist"))) stop("object must be of class 'pcrfit', 'modlist' or 'replist'!")
 if (CLASS[1] == "pcrfit") OBJ <- list(object) else OBJ <- object
 if (CLASS[2] == "replist") {
    OBJ <- rep2mod(object)
    GROUP <- attr(OBJ, "group")
 }
 
 outOBJ <- vector("list", length = length(OBJ))

 for (i in 1:length(OBJ)) {
     tempOBJ <- OBJ[[i]]
     NAMES <- tempOBJ$names
     if (verbose) cat(NAMES, "\nCalculating first and second derivative maximum...\n")
     flush.console()
     EFF <- try(efficiency(tempOBJ, plot = FALSE, ...), silent = TRUE)
     cpD2 <- if (inherits(EFF, "try-error")) 0 else EFF$cpD2
     cpD1 <- if (inherits(EFF, "try-error")) 0 else EFF$cpD1
     if (verbose) cat("Calculating R-square...\n")
     flush.console()
     RSQ <- if(inherits(EFF, "try-error")) 0 else Rsq(tempOBJ)
     if (verbose) cat("Checking for sigmoidal consistency...\n")

     if (cpD2 > cpD1 |  cpD1 - cpD2 > 10 | RSQ < 0.9) {
        tempOBJ$outlier <- TRUE
        if (verbose) cat(" Found non-sigmoidal structure for", NAMES, "...\n", sep = " ")
        flush.console()

        if (remove) {
           if (verbose) cat(" Removing", NAMES, "...\n\n", sep = " ")
           flush.console()
           next
        }
        
        if (verbose) cat(" Tagging name of", NAMES, "...\n", sep = " ")
        flush.console()
        tempOBJ$names <- paste("**", tempOBJ$names, "**", sep = "")
     } else tempOBJ$outlier <- FALSE
     
     cat("\n")
     outOBJ[[i]] <- tempOBJ
 }
 
 if (remove) {
  WHICH <- which(outOBJ == "NULL")
  if (length(WHICH) > 0) {
    outOBJ <- outOBJ[-WHICH]
    if (CLASS[2] == "replist") GROUP <- GROUP[-WHICH]
  }
 }      

 if (CLASS[1] == "pcrfit") outOBJ <- outOBJ[[1]]  
 
 class(outOBJ) <- CLASS 
 
 if (CLASS[2] == "replist") {    
    cat("Updating 'replist':\n")     
    outOBJ <- replist(outOBJ, GROUP, verbose = verbose, ...)
 }    

 return(outOBJ)
} 
