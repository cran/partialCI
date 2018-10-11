# partialCI
R package for fitting the partially cointegrated model

A collection of time series is partially cointegrated if a linear combination of these time series can be found so that the residual spread is partially autoregressive - meaning that it can be represented as a sum of an autoregressive series and a random walk. This concept is useful in modeling certain sets of financial time series and beyond, as it allows for the spread to contain transient and permanent components alike. Partial cointegration has been introduced by Clegg and Krauss (2017) \doi{10.1080/14697688.2017.1370122}, along with a large-scale empirical application to financial market data. The partialCI package comprises estimation, testing, and simulation routines for partial cointegration models in state space. Clegg et al. (2017) <http://hdl.handle.net/10419/150014> provide an in in-depth discussion of the package functionality as well as illustrating examples in the fields of finance and macroeconomics.
If a collection of time series is partially cointegrated, then the spread between
them can be interepreted as a mean-reverting series that has possibly been 
contaminated with a (hopefully small) random walk.

The developer version can be find on Github. To use the developer version of the partialCI package, you will need to start by installing it,
which can be done using devtools:

```
> install.packages("devtools")   # if devtools is not already installed
> install_github("matthewclegg/partialCI")
```

To find the partially cointegrated model that best fits two series 
X and Y, use:

```
> fit.pci(Y, X)
```

An interface to Yahoo! Finance permits you to find the best fits for
two particular stocks of interest:

```
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
```

This example was run on 29/9/2018.  RDS-A and RDS-B are two 
classes of shares offered by Royal Dutch Shell that differ slightly
in aspects of their tax treatment.  The above fit shows that
the spread between the two shares is mostly mean-reverting but that
it contains a small random walk component.  The mean-reverting
component accounts for 44.6% of the variance of the daily returns.
The value of 0.1770 for rho corresponds to a half-life of mean
reversion of about 4 trading days.

To test the goodness of fit, the test.pci function can be used:

```
> h <- yfit.pci("RDS-B", "RDS-A")
> test.pci(Y=h)

LR test of [RW or CI(1)] vs Almost PCI(1) (joint penalty,wilk)

data:  RDS-B / RDS-A

Hypothesis      Statistic  p-value    alpha alpha_bonf alpha_holm
Random Walk         -8.30    0.000    0.050    0.025    0.050
AR(1)               -7.68    0.000    0.050    0.025    0.025

```

The test.pci function tests each of two different null hypotheses:
(a) the residual series is purely a random walk, and (b) the residual series is
purely autoregressive. Only if both null hypothesis can be rejected a time series is classified as partially cointegrated. The two p-values of 0.000 indicate that RDS-A and RDS-B are indeed partially cointegrated. This still holds true if a Bonferroni corrected significance level (alpha_bonf) is considered.

The partialCI package also contains a function for searching for
hedging portfolios.  Given a particular stock (or time series),
a search can be conducted to find the set of stocks that best
replicate the target stock.  In the following example, a hedge 
is sought for SPY using sector ETF's.

```
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
```

The top table displays the quality of the fit that is found as each new
factor is added to the fit.  The best fit consisting of only one factor
is found by using XLK (the technology sector).  The negative log likelihod
score for this model is 1307.48.  However, the random walk
hypothesis (p[rw]) cannot be rejected at the 5% level.  When adding
XLI (the industrial sector), the negative log likelihood drops to 747.43
and the random walk and the purely autoregressive hypothesis can not be rejected at the 5% level if if we account for multiple testing. The best overall fit is obtained by also adding XLY (consumer discretionary) to the hedging
portfolio.  The final fit is

```
  SPY = 1.3450 XLK + 1.2436 XLI + 0.6750 XLY
```

For this fit, the proportion of variance attributable to the mean reverting
component is 41.0%. In addition, for the best fit we can reject the random walk and the purely autoregressive hypothesis at the 5% level. The later holds even if we account for multiple testing. 

Please feel free to write to us if you have questions or suggestions.

Jonas Rende 

jonas.rende@fau.de

Matthew Clegg  

matthewcleggphd@gmail.com  

Christopher Krauss

christopher.krauss@fau.de



Oktober 02, 2018  

