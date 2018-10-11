\name{partialCI-package}
\alias{partialCI-package}
\alias{partialCI}
\docType{package}
\title{
Partial Cointegration
}
\description{
A collection of time series is said to be partially cointegrated if they have a linear
combination that is partially autoregressive, e.g., that can be represented as a sum of an
autoregressive series and a random walk.  This may be useful in modeling
certain sets of financial time series.

To find the partially cointegrated model that best fits two series X and Y, use:

\preformatted{
> fit.pci(Y, X)
}

An interface to Yahoo! Finance permits you to find the best fits for two particular 
stocks of interest:

\preformatted{
> yfit.pci("RDS-B", "RDS-A")
Fitted values for PCI model
  Y[t] = X[t] %*% beta + M[t] + R[t]
  M[t] = rho * M[t-1] + eps_M [t], eps_M[t] ~ N(0, sigma_M^2)
  R[t] = R[t-1] + eps_R [t], eps_R[t] ~ N(0, sigma_R^2)

           Estimate Std. Err
beta_RDS-A   1.0427   0.0098
rho          0.1770   0.1715
sigma_M      0.0825   0.0130
sigma_R      0.1200   0.0095

-LL = -255.97, R^2[MR] = 0.446
}

This example was run on 29/9/2018.  RDS-A and RDS-B are two 
classes of shares offered by Royal Dutch Shell that differ slightly
in aspects of their tax treatment.  The above fit shows that
the spread between the two shares is mostly mean-reverting but that
it contains a small random walk component.  The mean-reverting
component accounts for 44.6\% of the variance of the daily returns.
The value of 0.1770 for rho corresponds to a half-life of mean
reversion of about 4 trading days.

To test the goodness of fit, the \code{test.pci} function can be used:

\preformatted{
> h <- yfit.pci("RDS-B", "RDS-A")
> test.pci(Y=h)

R test of [RW or CI(1)] vs Almost PCI(1) (joint penalty,wilk)

data:  RDS-B / RDS-A

Hypothesis      Statistic  p-value    alpha alpha_bonf alpha_holm
Random Walk         -8.30    0.000    0.050    0.025    0.050
AR(1)               -7.68    0.000    0.050    0.025    0.025
}

The \code{test.pci}  function tests each of two different null hypotheses:
(a) the residual series is purely a random walk, and (b) the residual series is
purely autoregressive. Only if both null hypothesis can be rejected a time series is classified as partially cointegrated. The two p-values of 0.000 indicate that RDS-A and RDS-B are indeed partially cointegrated. This still holds true if a Bonferroni corrected significance level (alpha_bonf) is considered.

The partialCI package also contains a function for searching for
hedging portfolios.  Given a particular stock (or time series),
a search can be conducted to find the set of stocks that best
replicate the target stock.  In the following example, a hedge 
is sought for SPY using sector ETF's.

The partialCI package also contains a function for searching for hedging portfolios. 
Given a particular stock (or time series), a search can be conducted to find the set 
of stocks that best replicate the target stock. In the following example, a hedge is 
sought for SPY using sector ETF's.

\preformatted{
 -LL   LR[rw]    p[rw]    p[mr]      rho  R^2[MR]   Factor |   Factor coefficients
 1307.48  -2.9362   0.0531   0.2847   0.9722   0.8362      XLK |   3.1505 
  747.43  -3.6282   0.0266   0.0253   0.9112   0.5960      XLI |   1.8554   1.6029 
  548.63  -5.9474   0.0026   0.0007   0.6764   0.4098      XLY |   1.3450   1.2436   0.6750 

Fitted values for PCI model
  Y[t] = X[t] %*% beta + M[t] + R[t]
  M[t] = rho * M[t-1] + eps_M [t], eps_M[t] ~ N(0, sigma_M^2)
  R[t] = R[t-1] + eps_R [t], eps_R[t] ~ N(0, sigma_R^2)

           Estimate Std. Err
beta_XLK     1.3450   0.0419
beta_XLI     1.2436   0.0352
beta_XLY     0.6750   0.0310
rho          0.6764   0.1502
sigma_M      0.2323   0.0429
sigma_R      0.3045   0.0331

-LL = 548.63, R^2[MR] = 0.410
}

The top table displays the quality of the fit that is found as each new
factor is added to the fit.  The best fit consisting of only one factor
is found by using XLK (the technology sector).  The negative log likelihod
score for this model is 1307.48.  However, the random walk
hypothesis (p[rw]) cannot be rejected at the 5\% level.  When adding
XLI (the industrial sector), the negative log likelihood drops to 747.43
and the random walk and the purely autoregressive hypothesis can not be rejected at the 5\% level if if we account for multiple testing. The best overall fit is obtained by also adding XLY (consumer discretionary) to the hedging
portfolio.  The final fit is

\preformatted{
  SPY = 1.3450 XLK + 1.2436 XLI + 0.6750 XLY
}

For this fit, the proportion of variance attributable to the mean reverting
component is 41.0\%. In addition, for the best fit we can reject the random walk and the purely autoregressive hypothesis at the 5\% level. The later holds even if we account for multiple testing. 

Please feel free to contact us if you have questions or suggestions.

Jonas Rende, Matthew Clegg and Christopher Krauss 

Oktober 02, 2018  
}
\author{
Jonas Rende \email{jonas.rende@fau.de}

Matthew Clegg \email{matthewcleggphd@gmail.com}

Christopher Krauss \email{christopher.krauss@fau.de}
}

\seealso{
\code{\link{fit.pci}} 
\code{\link{yfit.pci}}
\code{\link{test.pci}} 
\code{\link{hedge.pci}}
\code{\link{yhedge.pci}}
}

\keyword{ ts }
\keyword{ models }% __ONLY ONE__ keyword per line
