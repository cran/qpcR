\name{qpcR}
\alias{qpcR}

\docType{package}

\title{Modelling and analysis of real-time PCR data}

\description{
Model fitting, optimal model selection and calculation of various features that are essential
in the analysis of quantitative real-time polymerase chain reaction (qPCR).
}

\details{

\tabular{ll}{
Package: \tab qpcR\cr
Type: \tab Package\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Version: \tab 1.0-6\cr
Date: \tab 2008-01-28\cr
License: \tab GNU GPL (version 2 or any later version)\cr
Depends: \tab R (>= 2.0.0), drc\cr
URL: \tab http://www.r-project.org; http://www.bioassay.dk; http://www.dr-spiess.de \cr
Packaged: \tab Thu Jan 17 2008; spiess\cr
Built: \tab R 2.6.0; ; 2008-01-28; windows\cr
}

Index:
\preformatted{

package:
qpcR Modelling and analysis of real-time PCR data.

Functions:
\link{confband} Confidence bands for a sigmoidal fit.
\link{eff} The amplification efficiency curve of a fitted object.
\link{efficiency} Calculation of qPCR efficiency and several other important qPCR parameters.
\link{mchoice} Selection of the best model by nested f-tests.
\link{pcrplot} Plotting PCR data with fitted curves.
\link{repform} Formatting of PCR data with or without replicates for analysis.
\link{resVar} Residual variance of a sigmoidal fit.
\link{Rsq} R-square value of a sigmoidal fit.
\link{sliwin} Calculation of PCR efficiency by the window-of-linearity method.
\link{expfit} Calculation of PCR efficiency by fitting of an exponential model and minimization of residual variance.
\link{expfit2} Calculation of PCR efficiency by fitting of an exponential model and outlier cycle detection.
\link{pcrbatch} Batch calculation of qPCR efficiency and several other important qPCR parameters.
\link{RMSE} Root-mean-squared-error of a sigmoidal fit.
\link{calib} Calculation of qPCR efficiency by dilution curve analysis.
\link{pcrpred} Value prediction from a fitted sigmoidal model.
\link{pcrGOF} Summarize measures for the goodness-of-fit.
\link{expcomp} Comparison of several sigmodal models by RMSE within the exponential region.
\link{pcrimport} Simple qPCR data import function (i.e. from text files, clipboard).
\link{pcrfit} Simple wrapper function for fast model fitting. 
\link{AICc} Akaike's second-order Information Criterion for small sample numbers.
\link{evidence} Evidence ratio for model comparison with AIC or AICc.
\link{modlist} Create sigmoidal models from a dataframe and coerce them into a list.
\link{pcropt1} Optimization of efficiency and F0 estimation by combinatorial cycle elimination of plateau and ground phase cycles.

Datasets:
\link{reps} qPCR dilution experiment with replicates (Roche Lightcycler).
\link{reps2} Another qPCR dilution experiment with replicates (Roche Lightcycler).
\link{rutledge} qPCR data from Rutledge et al (Biorad Opticon).
\link{S27} qPCR data from 15 different testicular biopsies (Roche Lightcycler).
\link{reps3} A qPCR dilution experiment with replicates (Stratagene MXPro3000P).
}

}

\author{
Andrej-Nikolai Spiess <a.spiess@uke.uni-hamburg.de>, Christian
Ritz <ritz@kvl.dk>

Maintainer: Andrej-Nikolai Spiess <a.spiess@uke.uni-hamburg.de>
}

\references{
%~~ Literature or other references for background information ~~
}

\keyword{package}
\keyword{models}
\keyword{nonlinear}

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~ 
% Maybe there could be a link to 'drc' here?
}

\examples{

}
