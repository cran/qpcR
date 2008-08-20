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
Version: \tab 1.1-5\cr
License: \tab GNU GPL (version 2 or any later version)\cr
Depends: \tab R (>= 2.0.0), drc\cr
URL: \tab http://www.dr-spiess.de; http://www.r-project.org; http://www.bioassay.dk \cr
Packaged: \tab Tue Jun 15 2008; spiess\cr
Built: \tab R 2.7.1; ; 2008-08-20; windows\cr
}

package:
qpcR Modelling and analysis of real-time PCR data.

Functions:\cr
\link{confband} Confidence bands for a sigmoidal fit.\cr
\link{eff} The amplification efficiency curve of a fitted object.\cr
\link{efficiency} Calculation of qPCR efficiency and several other important qPCR parameters.\cr
\link{mchoice} Selection of the best model by nested F-tests, likelihood ratios or Akaike weights.\cr
\link{pcrplot} Plotting PCR data with fitted curves.\cr
\link{repform} Formatting of PCR data with or without replicates for analysis.\cr
\link{resVar} Residual variance of a fitted model.\cr
\link{Rsq} R-square value of a fitted model.\cr
\link{sliwin} Calculation of PCR efficiency by the window-of-linearity method.\cr
\link{expfit} Calculation of PCR efficiency by fitting an exponential model.\cr
\link{pcrbatch} Batch calculation of qPCR efficiency and several other important qPCR parameters.\cr
\link{RMSE} Root-mean-squared-error of a fitted model.\cr
\link{calib} Calculation of qPCR efficiency by dilution curve analysis.\cr
\link{pcrpred} Value prediction from a fitted sigmoidal model.\cr
\link{pcrGOF} Summarize measures for the goodness-of-fit.\cr
\link{expcomp} Comparison of several sigmodal models by RMSE within the exponential region.\cr
\link{pcrimport} Simple qPCR data import function (i.e. from text files, clipboard).\cr
\link{pcrfit} Simple wrapper function for fast model fitting.\cr
\link{AICc} Akaike's second-order corrected Information Criterion for small sample sizes.\cr
\link{evidence} Evidence ratio for model comparison with AIC or AICc.\cr
\link{modlist} Create nonlinear models from a dataframe and coerce them into a list.\cr
\link{pcropt1} Combinatorial elimination of plateau and ground phase cycles.\cr
\link{propagate} General error propagation function using covariates and Monte Carlo simulations.\cr
\link{ratiocalc} Calculation of ratios/propagated errors/t-statistics for all permutations/combinations of qPCR runs with or without reference data.\cr
\link{ratioplot} Plot method for data of class 'ratiocalc' which displays barplots of ratios with error bars.\cr
\link{curvemean} Building a model which averages a batch of qPCR curves.\cr
\link{PRESS} Allen's PRESS (Prediction Sum-Of-Squares) statistic.\cr
\link{akaike.weights} Calculation of Akaike weights/relative likelihoods/delta-AICs.\cr
\link{LR} Calculation of likelihood ratios for nested models.\cr
\link{RSS} Residual sum-of-squares of a fitted model.\cr
\link{residplot} A residuals bar-plot with colour-coded bars.\cr

Datasets:\cr
\link{reps} qPCR dilution experiment with replicates (Roche Lightcycler).\cr
\link{reps2} Another qPCR dilution experiment with replicates (Roche Lightcycler).\cr
\link{rutledge} qPCR data from Rutledge et al (Biorad Opticon).\cr
\link{S27} qPCR data from 15 different testicular biopsies (Roche Lightcycler).\cr
\link{reps3} A qPCR dilution experiment with replicates (Stratagene MXPro3000P).\cr
\link{guescini1} A qPCR dilution experiment with replicates (Lightcycler 480) from Guescini et al.\cr
\link{guescini2} A qPCR experiment with replicates and simulated inhibition (Lightcycler 480) from Guescini et al.\cr
}

\author{
Andrej-Nikolai Spiess <a.spiess@uke.uni-hamburg.de>, Christian
Ritz <ritz@kvl.dk>

Maintainer: Andrej-Nikolai Spiess <a.spiess@uke.uni-hamburg.de>
}

\keyword{package}
\keyword{models}
\keyword{nonlinear}



