.onAttach <- function(libname, pkgname)
{
    cat("\n")
    cat("'qpcR' has been loaded.\n\n")
    cat("Please cite R and the following if used for a publication:\n\n")
    cat("Spiess AN, Feig C, Ritz\nHighly accurate sigmoidal fitting of real-time PCR data by introducing a parameter for asymmetry.\nBMC Bioinformatics 2008, 29:221\n")
    cat("or\n")
    cat("Ritz C, Spiess AN. qpcR: an R package for sigmoidal model selection in quantitative real-time polymerase chain reaction analysis.\nBioinformatics 2008, 24:1549-1551\n\n") 
    cat("Newest version always available at www.dr-spiess.de/qpcR.html.\n\n")
}
