ABImport <- function(path = NULL, skip = 0, well.pos = 1, dye.pos = 2, ref.dye = "ROX")
{
  ## if path is not given, open selection window
  ## define path with double slashes, i.e. c:\\temp\\name.csv
  if (is.null(path)) cf <- choose.files() else cf <- path

  ## read in data, maybe skipping lines
  DATA <- read.csv(cf, skip = skip)  
  
  ## get unique well numbers 
  well.num <- unique(DATA[, well.pos])
   
  ## get unique dye names
  dye.names <- unique(DATA[, dye.pos])   
  
  ## check if reference dye exists
  ## and use other names as dye names
  hasREF <- ref.dye %in% dye.names   
  if(hasREF) {   
    m <- match(ref.dye, dye.names)
    dye.names <- dye.names[-m] 
  }
  
  ## get columns with fluorescence data
  fluo.col <- (dye.pos + 1):ncol(DATA)   
  
  ## create lists and dataframes to store the processed data
  origLIST <- NULL
  normLIST <- NULL
  origFRAME <- NULL
  normFRAME <- NULL 
    
  ## Normalize (if ref dye is present) non-reference dyes and put 
  ## into new data frame along with the original data     
  for (i in 1:length(dye.names)) {      
    for (j in well.num) {
      DYE <- which(DATA[, well.pos] == j & DATA[, dye.pos] == dye.names[i])
      if (hasREF) REF <- which(DATA[, well.pos] == j & DATA[, dye.pos] == ref.dye)
      origDAT <- DATA[DYE, fluo.col]        
      if (hasREF) normDAT <- origDAT/DATA[REF, fluo.col] 
      origFRAME <- cbind(origFRAME, t(origDAT))         
      if (hasREF) normFRAME <- cbind(normFRAME, t(normDAT))                          
    }          
    ## create column names from original well data
    if (hasREF) colnames(normFRAME) <- paste("Well.", well.num, sep = "")
    colnames(origFRAME) <- paste("Well.", well.num, sep = "") 
    ## add Cycle numbers in first column
    if (hasREF) normFRAME <- cbind(Cycles = 1:nrow(normFRAME), normFRAME)
    origFRAME <- cbind(Cycles = 1:nrow(origFRAME), origFRAME)       
    ## put into list
    if (hasREF) normLIST[[i]] <- normFRAME
    origLIST[[i]] <- origFRAME
    ## clear dataframe
    origFRAME <- NULL  
    normFRAME <- NULL   
  }
  ## give list item dye names
  names(origLIST) <- dye.names 
  if (hasREF) names(normLIST) <- dye.names 
  ## return lists from function
  if (!hasREF) normLIST <- origLIST 
  return(list(orig = origLIST, norm = normLIST))   
}   