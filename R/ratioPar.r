ratioPar <- function(
group = NULL, 
effVec = NULL,
cpVec = NULL,
type.eff = "individual",
plot = TRUE,
combs = c("same", "across", "all"),
refmean = FALSE,
verbose = TRUE,
...)
{
    ## create dummy "pcrbatch" data
    DAT <- matrix(nrow = 40, ncol = length(group) + 1)
    colnames(DAT) <- paste("Run", 1:ncol(DAT), sep = ".")
    colnames(DAT)[1] <- "Vars"
    DAT <- data.frame(DAT)
    class(DAT)[2] <- "pcrbatch"
    
    ratiobatch(data = DAT, group = group, which.eff = effVec, which.cp = cpVec,
               type.eff = type.eff, plot = plot, combs = combs, refmean = refmean,
               verbose = verbose, ...)
}