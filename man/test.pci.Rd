\name{test.pci}

\alias{test.pci}
\title{Tests the goodness of fit of a partial cointegration model}
\description{
Tests the goodness of fit of a partial cointegration model
}
\usage{
test.pci(Y, X,
        pci_opt_method=c("jp", "twostep"),
        irobust=FALSE,
        inu  = 5,
        null_hyp=c("par","rw", "ar1"),
        imethod = c("wilk","boot"),
        inrep = 999,
        istart.seed= 1,
        alpha = 0.05,
        use.multicore = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  The time series that is to be modeled.  A plain or \code{\link{zoo}} vector of length \code{n} or a \code{\link{fit.pci}} object.
}
  \item{X}{
  A (possibly \code{\link{zoo}}) matrix of dimensions \code{n} x \code{k}.  
  If \code{k=1}, then this may be a plain or \code{\link{zoo}} vector.
}

\item{pci_opt_method}{
The method that will be used for fitting a partially cointegrated
model to \code{X} and \code{Y}. This can be either \code{"jp"}
(joint penalty) or \code{"twostep"} (Engle-Granger two-step).
See \code{\link{fit.pci}} for a complete explanation.
Default: \code{"jp"}.
}

\item{irobust}{
  If \code{TRUE}, then the residuals are assumed to follow a t-distribution.  Default: \code{FALSE}.
}

 \item{inu}{
The degrees-of-freedom parameter to be used in robust estimation.  Default: 5.
}
  
  \item{null_hyp}{
This specifies the null hypothesis.
This can be either \code{"rw"}, \code{"ar1"} or \code{"par"}.
If \code{"rw"}, then the null hypothesis is a random walk.
If \code{"ar1"}, then the null hypothesis is an autoregressive process
of order 1.  (In this case, the null hypothesis calls for Y and X to be
cointegrated.)  If \code{"par"}, then the null hypothesis is that the process can be modelled as a PAR process.  Default: "par".
}

\item{imethod}{
The method used to calculate p-values associated with the likelihood ratio test. If set to "wilk" p-values are determined under the assumption that theorem of Wilks holds (Wilks, 1933). If set to "boot" a parametric bootstrap method is carried out following the instructions of MacKinnon (2009). For details see description below.
Default: \code{"wilk"}.
}

\item{inrep}{
The number of bootstrap iterations. Only han an effect if code{imethod = c("boot")}. Default: inrep = 999.
}

\item{istart.seed}{
Sets the starting seed for the parametric bootstrap. Only han an effect if code{imethod = c("boot")}. Default: istart.seed= 1.
}



\item{alpha}{
The global significance level used to determine the local significance level following (i) the conservative Bonferroni correction and (ii) the more liberal Holm correction to account for multiple testing. Only important if \code{null_hyp="par"}. Default: 0.05.
}

\item{use.multicore}{
If \code{TRUE}, parallel processing will be used to improve performance. Only han an effect if  code{imethod = c("boot")}. Default: \code{FALSE}.
}

}
\details{
The likelihood ratio test is used to determine whether 
the null hypothesis should be rejected or not, since the null models are nested in the partial cointegration alternative (Neymann, 1933).  That is to say, a search is performed for the best fitting model under the null hypothesis, and the log likelihood score of this model is computed.
Then a search is performed for the best fitting model under the alternative hypothesis of partial cointegration, and the log
likelihood score of this model is computed. We consider two null hypotheses, namely the the pure random walk hypotheses and the pure stationary AR(1) hypothesis (refers to classic cointegration). Given that under the null hypothesis we are solely testing parameters at the boundaries of the corresponding parameter space the regularity condition underlying the theorem of Wilks that the true parameter has to be an inner point of the parameter space does not hold (Wilks, 1938). At boundaries standard asymptotics fail. The asymptotic distribution of the likelihood ratio test statistic is a mixture of \eqn{\chi^2} distributions (Shapiro, 1988): 
\deqn{w[0] * \chi[0]^2 + w[1] * \chi[1]^2 + ... + w[m] * \chi[m]^2 + ... +  w[k] * \chi[k]^2}
where the weights \eqn{w[i]}, with \eqn{i = \{0, ... ,k\}} sum up to one. Let m denote the degrees of freedom associated with the \eqn{\chi^2} distribution of interest and \eqn{k} is equal to the parameter difference of the full and the null model. Note that \eqn{\chi^2[0] = 0}, indicating a point mass at zero. Under Wilks theorem the weight for the \eqn{\chi^2[k]} distribution is set to one. Therefore, the resulting p-values falsely assuming Wilks theorem holds are more conservative compared to the p-values associated with a mixed distribution. Therefore, we implement the latter as a conservative approximation (Stoel, 2006) for the true underlying distribution \code{(imethod = "wilk")}. If \code{(imethod = "boot")}, a parametric bootstrap following the instructions of MacKinnon (2009) is carried out the calculate the p-value of the corresponding null hypothesis.  
}
\value{An object of class \code{"pcitest"} containing the
results of the hypothesis test.

}

\references{
Clegg, Matthew, 2015.
  Modeling Time Series with Both Permanent and Transient Components
  using the Partially Autoregressive Model.
\emph{Available at SSRN: http://ssrn.com/abstract=2556957}

 Clegg, Matthew and Krauss, Christopher, 2018. Pairs trading with partial cointegration. Quantitative Finance, 18(1). 
 
 MacKinnon, J. G., 2009. Bootstrap Hypothesis Testing: 6. In Handbook of Computational Econometrics, Wiley-Blackwell
 
  Neyman, J. and Pearson, E. S.,1933. On the Problem of the Most Efficient Tests of Statistical Hypotheses. Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character,231.
  
    Shapiro, A., 1988. Towards a Unified Theory of Inequality Constrained Testing in Multivariate Analysis. International Statistical Review / Revue Internationale de Statistique,56(1):49
     Stoel, R. D., Garre, F. G., Dolan, C., and van den Wittenboer, G., 2006. On the likelihood ratio test in structural equation modeling when parameters are subject to boundary
constraints. Psychological methods, 11(4).
   
  

}

\author{
Matthew Clegg \email{matthewcleggphd@gmail.com}

Christopher Krauss \email{christopher.krauss@fau.de}

Jonas Rende \email{jonas.rende@fau.de}
}
\seealso{
\code{\link{fit.pci}} Fits a partially cointegrated model

\code{\link{likelihood_ratio.pci}} Computes the likelihood ratio of a PCI model
  versus a null model
}
\examples{
# The following should reject both the random walk and AR(1) models

\dontrun{
# Example using the very fast "wilk" method
# The following should be classified as PCI, i.e., both Nullhypothesis should be rejected even
# if we account for multiple testing (alpha_bonf, alpha_holm)
set.seed(1313)
YX <- rpci(n=1000,beta=c(1), sigma_C=c(0.1), rho=0.8, sigma_M=1, sigma_R=1)
test.pci(Y=YX[,1], X=YX[,2:ncol(YX)])

# Example using the parametric bootstrap "boot" method with 999 iterations (inrep = 999 and 
# a starting seed of 1 (istart.seed= 1)
# The following should be classified as PCI, i.e., both Nullhypothesis should be rejected even
# if we account for multiple testing (alpha_bonf, alpha_holm)
# Results are very similar
set.seed(1313)

YX <- rpci(n=1000,beta=c(1), sigma_C=c(0.1), rho=0.8, sigma_M=1, sigma_R=1)
test.pci(Y=YX[,1], X=YX[,2:ncol(YX)],imethod = "boot", inrep = 999,
        istart.seed= 1)
        

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ models }
