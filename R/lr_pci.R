
# Functions for computing likelihoods, likelihood ratios and likelihood ratio
# tests for partially cointegrated (PCI) series

# Copyright (C) 2016 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

utils::globalVariables(c("PCI.SAMPLES", "PCI.SAMPLES.DT"))

pci.nu.default <- function () 5

par.nu.default <- function () {
  # The default value to be used for the degrees of freedom parameter
  # when using robust estimation.
  
  5
}




loglik.par.kfas <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1]) {
  # Given a sequence Y and a parameterization (rho, sigma_M, sigma_R) of an
  # associated PAR process, calculates the negative log likelihood that
  # Y would be observed under these process parameters.  
  
  if (length(Y) < 1) return(NA_real_)
  if (length(dim(Y)) > 0) Y <- Y[,1]
  Y <- zoo::coredata(Y)
  
  Zt <- matrix(c(1, 1), 1, 2)
  Ht <- matrix(c(0), 1, 1)
  Tt <- matrix(c(rho, 0, 0, 1), 2, 2, byrow = TRUE)
  Rt <- matrix(c(1, 0, 0, 1), ncol=2, nrow=2, byrow = TRUE)
  a1 <- matrix(c(M0, R0), 2, 1)
  Qt <- matrix(c(sigma_M^2, 0, 0 , sigma_R^2), 2, 2)
  P1 <- matrix(c(sigma_M^2, 0, 0, sigma_R^2), 2, 2,byrow = TRUE)
  P1inf<-matrix(c(0, 0, 0, 0), 2, 2,byrow = TRUE)
  model_par<-KFAS::SSModel(Y ~ -1 + KFAS::SSMcustom(Z = Zt, T = Tt, Q = Qt, a1 = a1, P1 = P1,
                                                    P1inf = P1inf), H = Ht)
  
  res<- (-1)*logLik(model_par, check.model = FALSE, convtol = 1e-08)
  return(res)
  
}

loglik.par.ss <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1]) {
  # Given a sequence Y and a parametrization (rho, sigma_M, sigma_R) of an 
  # associated PAR process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter.
  #
  # Despite the fact that this is hand-coded, the likelihood function that uses
  # the fast kalman filter implementation given above is about twice as fast,
  # even after turning on the byte compiler.  
  
  if (length(dim(Y)) > 0) Y <- Y[,1]
  Y <- coredata(Y)
  
  M0 <- as.numeric(M0)
  R0 <- as.numeric(R0)
  
  n <- length(Y)
  if (n < 1) return(NA_real_)
  
  K <- kalman.gain.par(rho, sigma_M, sigma_R)
  esumsq <- 0
  M <- M0
  R <- R0
  tvar <- sigma_M^2 + sigma_R^2
  
  for (i in 1:n) {
    xhat <- rho * M + R
    e <- Y[i] - xhat
    esumsq <- esumsq + e*e
    M <- rho * M + e * K[1]
    R <- R + e * K[2]
  }
  
  nll <- (n/2)*log(tvar * 2*pi) + esumsq/(2*tvar)
  
  nll
}

llst <- function(X, mu, sigma, df) {
  # Computes the negative log likelihood that (X - mu)/sigma is distributed
  # according to a t-distribution with df degrees of freedom.
  # Input values:
  #   X:     Vector of observations whose likelihood is to be computed
  #   mu:    The location parameter of the distribution
  #   sigma: The scale parameter -- similar to standard deviation
  #   df:    The degrees of freedom of the distribution
  #
  # This likelihood function has been optimized by computing the
  # values of the gamma function and some other constants only once,
  # rather than re-computing these values for each element of X.
  # This results in a factor of 10 (or more) speedup over the
  # naive implementation.
  
  if ((sigma <= 0) || (df <= 1) || (length(X) < 1)) return (NA_real_)
  
  # Following is the non-optimized version of the code:
  #  D <- dt((X - mu) / sigma, df) / sigma
  #  ll <- -sum(log(D))
  #  return(ll)
  
  # Following is an equivalent optimized version:
  N <- length(X)
  C <- N * (lgamma((df+1)/2) - 0.5*log(df * pi) - lgamma(df/2) - log(sigma))
  S <- -((df + 1)/2) * sum(log(1 + ((X-mu)/sigma)^2/df))
  ll <- -(S+C)
  ll
}

loglik.par.ss.t <- function (Y, rho, sigma_M, sigma_R, M0=0, R0=Y[1], nu=par.nu.default()) {
  # Given a sequence Y and a parametrization (rho, sigma_M, sigma_R) of an 
  # associated PAR process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter.
  
  if (length(dim(Y)) > 0) Y <- Y[,1]
  Y <- coredata(Y)
  
  M0 <- as.numeric(M0)
  R0 <- as.numeric(R0)
  
  K <- kalman.gain.par(rho, sigma_M, sigma_R)
  if (is.na(K[1])) return(NA_real_)
  
  n <- length(Y)
  if (n < 1) return(NA_real_)
  esumsq <- 0
  M <- M0
  R <- R0
  tvar <- sigma_M^2 + sigma_R^2
  tsd <- sqrt(tvar)
  
  evec <- numeric(n)
  for (i in 1:n) {
    xhat <- rho * M + R
    e <- Y[i] - xhat
    evec[i] <- e
    M <- rho * M + e * K[1]
    R <- R + e * K[2]
  }
  
  #    nll <- llst(evec, 0, tsd * sqrt((nu - 2) / nu), nu)
  nll <- llst(evec, 0, tsd, nu)
  
  nll
}








loglik.pci.kfas<- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0) {
  # Given a sequence Y, a basis X, and a parameterization 
  # (beta, rho, sigma_M, sigma_R) of an
  # associated PCI process, calculates the negative log likelihood that
  # Y would be observed under these process parameters.  
  
  if (is.null(dim(X)) && (length(beta) == 1)) {
    Z <- as.numeric(Y - X * beta)
  } else {
    Z <- as.numeric(Y - X %*% beta )
  }
  if (missing(R0)) R0 <- Z[1]
  loglik.par.kfas (Z, rho, sigma_M, sigma_R, M0, R0)
}

loglik.pci.ss <- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0) {
  # Given a sequence Y, basis X, and a parametrization 
  # (beta, rho, sigma_M, sigma_R) of an 
  # associated PCI process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter.
  
  if (is.null(dim(X)) && (length(beta) == 1)) {
    Z <- as.numeric(Y - X * beta)
  } else {
    Z <- as.numeric(Y - X %*% beta)
  }
  if (missing(R0)) R0 <- Z[1]
  loglik.par.ss (Z, rho, sigma_M, sigma_R, M0, R0)
}

loglik.pci.css <- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0) {
  # Given a sequence Y, basis X, and a parametrization 
  # (beta, rho, sigma_M, sigma_R) of an 
  # associated PCI process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter.  Uses a C implementation.
  
  if (is.null(dim(X)) && (length(beta) == 1)) {
    Z <- as.numeric(Y - X * beta)
  } else {
    Z <- as.numeric(Y - X %*% beta)
  }
  if (missing(R0)) R0 <- Z[1]
  loglik_par_c (Z, rho, sigma_M, sigma_R, M0, R0)
}

loglik.pci.sst <- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0, nu=pci.nu.default()) {
  # Given a sequence Y, basis X, and a parametrization 
  # (beta, rho, sigma_M, sigma_R) of an 
  # associated PCI process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter and a robust probability distribution.
  
  if (is.null(dim(X)) && (length(beta) == 1)) {
    Z <- as.numeric(Y - X * beta)
  } else {
    Z <- as.numeric(Y - X %*% beta)
  }
  if (missing(R0)) R0 <- Z[1]
  loglik.par.ss.t (Z, rho, sigma_M, sigma_R, M0, R0, nu=nu)
}

loglik.pci.csst <- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0, nu=pci.nu.default()) {
  # Given a sequence Y, basis X, and a parametrization 
  # (beta, rho, sigma_M, sigma_R) of an 
  # associated PCI process, calculates the negative log likelihood that Y would
  # be observed under these process parameters, using a steady state 
  # Kalman filter and a robust probability distribution.  Uses a C implementation.
  
  if (is.null(dim(X)) && (length(beta) == 1)) {
    Z <- as.numeric(Y - X * beta )
  } else {
    Z <- as.numeric(Y - X %*% beta)
  }
  if (missing(R0)) R0 <- Z[1]
  loglik_par_t_c (Z, rho, sigma_M, sigma_R, M0, R0, nu)
}

loglik.pci <- function (Y, X, beta, rho, sigma_M, sigma_R, M0=0, R0=0, 
                        calc_method=c("css", "kfas", "ss", "sst", "csst"), nu=pci.nu.default()) {
  # Given a sequence Y, basis X, and a parameterization (beta, rho, sigma_M, sigma_R) of an
  # associated PCI process, calculates the negative log likelihood that
  # Y would be observed under these process parameters.  The method used
  # for calculating the log likelihood is determined by "par_method":
  #   kfas:  Uses the Fast Kalman Filter (KFAS) package
  #   ss:   Uses a steady state Kalman filter
  #   css:  Uses a steady state Kalman filter coded in C
  
  switch(match.arg(calc_method),
         kfas=loglik.pci.kfas(Y, X, beta, rho, sigma_M, sigma_R, M0, R0),
         ss=loglik.pci.ss(Y, X, beta, rho, sigma_M, sigma_R, M0, R0),
         css=loglik.pci.css(Y, X, beta, rho, sigma_M, sigma_R, M0, R0),
         sst=loglik.pci.sst(Y, X, beta, rho, sigma_M, sigma_R, M0, R0, nu=nu),
         csst=loglik.pci.csst(Y, X, beta, rho, sigma_M, sigma_R, M0, R0)
  )
  
}

likelihood_ratio.pci <- function (
  Y,                       # The series which is being fit
  X,                       # The basis used for hedging X 
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  null_model=c("rw", "ar1"),  # Specifies the null hypothesis
  # rw = null model estimates sigma_R, assuming rho = sigma_M = 0.
  #      This is the default.
  # ar1 = null model estimates rho and sigma_M, assuming sigma_R = 0.
  pci_opt_method=c("jp", "twostep"),
  # Method to be used for fitting Y to X.
  #   jp:  The coefficients of Y are jointly optimized
  #     with the parameters of the AAR fit of the residuals
  #   twostep: A modified Engle-Granger procedure is used, where
  #     the coefficients of Y are first estimated, and then an AAR
  #     model is fit to the residuals.
  nu=5                     # If robust is TRUE, the degrees of freedom parameter                          
) {
  null_model <- match.arg(null_model)
  pci_opt_method <- match.arg(pci_opt_method)
  
  f.alt <- fit.pci(Y, X, robust=robust, pci_opt_method=pci_opt_method, nu=nu)
  f.null <- fit.pci(Y, X, robust=robust, par_model=null_model, pci_opt_method=pci_opt_method, nu=nu)
  f.alt$negloglik - f.null$negloglik   
}

likelihood_ratio.pci_internal <- function (
  Y,                       # The series which is being fit
  X,                       # The basis used for hedging X 
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  null_model=c("rw", "ar1"),  # Specifies the null hypothesis
  # rw = null model estimates sigma_R, assuming rho = sigma_M = 0.
  #      This is the default.
  # ar1 = null model estimates rho and sigma_M, assuming sigma_R = 0.
  pci_opt_method=c("jp", "twostep"),
  # Method to be used for fitting Y to X.
  #   jp:  The coefficients of Y are jointly optimized
  #     with the parameters of the AAR fit of the residuals
  #   twostep: A modified Engle-Granger procedure is used, where
  #     the coefficients of Y are first estimated, and then an AAR
  #     model is fit to the residuals.
  nu=5                     # If robust is TRUE, the degrees of freedom parameter                          
) {
  null_model <- match.arg(null_model)
  pci_opt_method <- match.arg(pci_opt_method)
  
  f.alt <- fit.pci(Y, X, robust=robust, pci_opt_method=pci_opt_method, nu=nu)
  f.null <- fit.pci(Y, X, robust=robust, par_model=null_model, pci_opt_method=pci_opt_method, nu=nu)
  lr<- f.alt$negloglik - f.null$negloglik  
  suppressWarnings(
   data.frame(likelihood_ratio = lr, 
             beta=f.alt$beta,
             rho=f.alt$rho,
             sigma_M=f.alt$sigma_M,
             sigma_R=f.alt$sigma_R,
             pvmr=f.alt$pvmr,
             beta_n=f.null$beta,
             rho_n=f.null$rho,
             sigma_M_n=f.null$sigma_M,
             sigma_R_n=f.null$sigma_R,
             pvmr_n=f.null$pvmr)
  )
}






print.pcitest <- function (x, ...) {

  print.pcitest.internal (x, ...)
}



print.pcitest.internal <- function (AT, ...) {
  # See stats:::print.htest
  cat("\n")
  cat(strwrap(AT$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", AT$data.name, "\n", sep = "")
  cat("\n")
  
  if (AT$robust) {
    h0a <- "Robust RW"
    h0b <- "Robust AR(1)"
    h1 <- "Robust PAR"
  } else {
    h0a <- "Random Walk"
    h0b <- "AR(1)"
    h1 <- "PAR"
  }
  
 
  cat(sprintf("%-12s %12s %8s %8s %8s %8s\n", "Hypothesis", "Statistic", "p-value","alpha","alpha_bonf", "alpha_holm"))
  if(AT$p.value[1]<AT$p.value[2]){
  cat(sprintf("%-12s %12.2f %8.4f %8.4f %8.3f %8.3f\n", h0a, AT$statistic[1], AT$p.value[1],AT$alpha,AT$alpha/2,AT$alpha/2))
  cat(sprintf("%-12s %12.2f %8.4f %8.4f %8.3f %8.3f\n", h0b, AT$statistic[2], AT$p.value[2],AT$alpha,AT$alpha/2,AT$alpha))
  } else {
    cat(sprintf("%-12s %12.2f %8.3f %8.3f %8.3f %8.3f\n", h0a, AT$statistic[1], AT$p.value[1],AT$alpha,AT$alpha/2,AT$alpha))
    cat(sprintf("%-12s %12.2f %8.3f %8.3f %8.3f %8.3f\n", h0b, AT$statistic[2], AT$p.value[2],AT$alpha,AT$alpha/2,AT$alpha/2))
  }
  
  cat("\n")
  invisible(AT)
}



test.pci<-function(Y,
                   X,
                   pci_opt_method=c("jp", "twostep"),
                   irobust=FALSE,
                   inu  = 5,
                   null_hyp=c("par","rw", "ar1"),
                   imethod = c("wilk","boot"),
                   inrep = 999,
                   istart.seed= 1,
                   alpha = 0.05,
                   use.multicore = F
){
  
  
  if (is(Y, "pci.hedge")) {
    Yhedge <- Y
    Y <- Yhedge$pci
  }
  
  if (is(Y, "pci.fit")) {
    Yfit <- Y
    Y <- Yfit$data
    X <- Yfit$basis
    if (missing(irobust)) irobust <- Yfit$robust
    if (missing(pci_opt_method)) {
      if(Yfit$pci.fit=="jointpenalty"){
      pci_opt_method <- "jp"
      } else {
        pci_opt_method <- "twostep"
      }
    }
  } 
  

  

  imethod <- match.arg(imethod)
  pci_opt_method <-  match.arg(pci_opt_method) 
  null_hyp <-match.arg(null_hyp)

  if(null_hyp=="rw"){
  p.rw<-test.pcarma.rw(X=X,Y=Y,method=imethod,nboot = inrep,iseed=istart.seed,irobust=irobust,nu  = inu,pci_opt_method =  pci_opt_method,use.multicore=use.multicore)
    return(p.rw)
  } else if(null_hyp=="ar1"){  
  p.ar<-test.pcarma.ar(X=X,Y=Y,method=imethod,nboot = inrep,iseed=istart.seed,irobust=irobust,nu  = inu,pci_opt_method =  pci_opt_method,use.multicore=use.multicore)
   return(p.ar)
   } else {
    p.rw<-test.pcarma.rw(X=X,Y=Y,method=imethod,nboot = inrep,iseed=istart.seed,irobust=irobust,nu  = inu,pci_opt_method =  pci_opt_method,use.multicore=use.multicore)
    p.ar<-test.pcarma.ar(X=X,Y=Y,method=imethod,nboot = inrep,iseed=istart.seed,irobust=irobust,nu  = inu,pci_opt_method =  pci_opt_method,use.multicore=use.multicore)
    
  }

  if (!irobust) {
    METHOD <- "LR test of [RW or CI(1)] vs Almost PCI(1)"
    alternative <- "Almost PCI(1)"
  } else {
    METHOD <- "LR test of [RW or Robust CI(1)] vs Robust Almost PCI(1)"
    alternative <- "Robust Almost PCI(1)"
  }
  
  if (pci_opt_method == "jp") {
    METHOD <- paste(METHOD, paste("(joint penalty,",imethod,")",sep=""))
  } else {
    METHOD <- paste(METHOD, paste("(two step,",imethod,")",sep=""))
  }
  STAT <- c(p.rw$statistic, p.ar$statistic)
  PVAL <- c(p.rw$p.value,p.ar$p.value)
  
  names(STAT) <- "LL"
  if(is.null(colnames(Y)) && !is.null(colnames(X))){
  name_buffer <- c()
  for(j in 1:length(colnames(X))){
  if(j==1){
  name_buffer[1] <- paste("Y", "/", colnames(X)[1])
  } else {
    name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
  } 
  }
  DNAME <- name_buffer[length(name_buffer)]
  } else if(!is.null(colnames(Y)) && is.null(colnames(X))){
    DNAME<-paste(colnames(Y),  "/", "X")
    
  } else if(!is.null(colnames(Y)) && !is.null(colnames(X))){
    name_buffer <- c()
    for(j in 1:length(colnames(X))){
      if(j==1){
        name_buffer[1] <- paste(colnames(Y), "/", colnames(X)[1])
      } else {
        name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
      } 
    }
    DNAME <- name_buffer[length(name_buffer)]
  } else {
    DNAME<-paste("Y",  "/", "X")
  }    
  structure(list(statistic = STAT, alternative = alternative, 
                 p.value = PVAL, method = METHOD, data.name = DNAME, robust=irobust,alpha=alpha),
            class = c("pcitest", "htest"))
  
 

}


test.pcarma.rw<-function(X,
                         Y,
                         method = c("wilk","boot"),
                         nboot = 999,
                         iseed = 1,
                         irobust=FALSE,
                         nu  = 5,
                         pci_opt_method="jp",
                         use.multicore=F
                              ){
  
  
  

  method <- match.arg(method)
  origin_fit<-likelihood_ratio.pci_internal(Y=Y,X=X,null_model = "rw",pci_opt_method=pci_opt_method, robust=irobust, nu=nu)
  if(method == "wilk"){
    p_value <- 1 - pchisq(-2 * origin_fit$likelihood_ratio[1], 2)

  } else if(method == "boot"){
    if(use.multicore){
      lrs <- c(parallel::mclapply(1:nboot,generate_null_model,lr_df = origin_fit,X_orig = X,inull_model ="rw",iseed=iseed,irobust=irobust,nu=nu,pci_opt_method=pci_opt_method),
               recursive=TRUE)
    } else {
      lrs <- c(lapply(1:nboot,generate_null_model,lr_df = origin_fit,X_orig = X,inull_model ="rw",iseed=iseed,irobust=irobust,nu=nu,pci_opt_method=pci_opt_method),
               recursive=TRUE)
    }
    lrs_sort <- sort(lrs)
    p_value<-length(which(lrs_sort<=origin_fit$likelihood_ratio[1]))/nboot

  }
  
  
 
  
  STAT <-   origin_fit$likelihood_ratio[1]
  PVAL <- p_value
  if (irobust) {
    METHOD <- "LR test of Robust RW vs Robust PCI(1)"
    alternative <- "RPCI(1)"
  } else {
    METHOD <- "LR test of RW vs PCI(1)"
    alternative <- "PCI(1)"
  }
  if (pci_opt_method == "jp") {
    METHOD <- paste(METHOD, paste("(joint penalty,",method,")",sep=""))
  } else {
    METHOD <- paste(METHOD, paste("(two step",method,")",sep=""))
  }
  names(STAT) <- "LL"
  if(is.null(colnames(Y)) && !is.null(colnames(X))){
    name_buffer <- c()
    for(j in 1:length(colnames(X))){
      if(j==1){
        name_buffer[1] <- paste("Y", "/", colnames(X)[1])
      } else {
        name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
      } 
    }
    DNAME <- name_buffer[length(name_buffer)]
  } else if(!is.null(colnames(Y)) && is.null(colnames(X))){
    DNAME<-paste(colnames(Y),  "/", "X")
    
  } else if(!is.null(colnames(Y)) && !is.null(colnames(X))){
    name_buffer <- c()
    for(j in 1:length(colnames(X))){
      if(j==1){
        name_buffer[1] <- paste(colnames(Y), "/", colnames(X)[1])
      } else {
        name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
      } 
    }
    DNAME <- name_buffer[length(name_buffer)]
  } else {
    DNAME<-paste("Y",  "/", "X")
  }     
  structure(list(statistic = STAT, alternative = alternative, 
                 p.value = PVAL, method = METHOD, data.name = DNAME),
            class = "htest")
  
  
 
}




test.pcarma.ar<-function(X,
                         Y,
                         method = c("wilk","boot"),
                         nboot = 999,
                         iseed = 1,
                         irobust=FALSE,
                         nu  = 5,
                         pci_opt_method="jp",
                         use.multicore=F
){

 
  method <- match.arg(method)
  origin_fit<-likelihood_ratio.pci_internal(Y=Y,X=X,null_model ="ar1",pci_opt_method=pci_opt_method, robust=irobust, nu=nu)
  if(method == "wilk"){
    p_value <- 1 - pchisq(-2 * origin_fit$likelihood_ratio[1], 1)

  } else if(method == "boot"){
 
  
    if(use.multicore){
      lrs <- c(parallel::mclapply(1:nboot,generate_null_model,lr_df = origin_fit,X_orig = X,inull_model ="ar1",iseed=iseed,irobust=irobust,nu=nu,pci_opt_method=pci_opt_method),
               recursive=TRUE)
    } else {

      lrs <- c(lapply(1:nboot,generate_null_model,lr_df = origin_fit,X_orig = X,inull_model ="ar1",iseed=iseed,irobust=irobust,nu=nu,pci_opt_method=pci_opt_method),
               recursive=TRUE)
    }
      lrs_sort <- sort(lrs)
      p_value<-length(which(lrs_sort<=origin_fit$likelihood_ratio[1]))/nboot
    }
    
  
  STAT <-   origin_fit$likelihood_ratio[1]
  PVAL <- p_value
  if (irobust) {
    METHOD <- "LR test of Robust CI vs Robust PCI(1)"
    alternative <- "RPCI(1)"
  } else {
    METHOD <- "LR test of CI vs PCI"
    alternative <- "PCI(1)"
  }
  if (pci_opt_method == "jp") {
    METHOD <- paste(METHOD, paste("(joint penalty,",method,")",sep=""))
  } else {
    METHOD <- paste(METHOD, paste("(two step",method,")",sep=""))
  }
  names(STAT) <- "LL"
  if(is.null(colnames(Y)) && !is.null(colnames(X))){
    name_buffer <- c()
    for(j in 1:length(colnames(X))){
      if(j==1){
        name_buffer[1] <- paste("Y", "/", colnames(X)[1])
      } else {
        name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
      } 
    }
    DNAME <- name_buffer[length(name_buffer)]
  } else if(!is.null(colnames(Y)) && is.null(colnames(X))){
    DNAME<-paste(colnames(Y),  "/", "X")
    
  } else if(!is.null(colnames(Y)) && !is.null(colnames(X))){
    name_buffer <- c()
    for(j in 1:length(colnames(X))){
      if(j==1){
        name_buffer[1] <- paste(colnames(Y), "/", colnames(X)[1])
      } else {
        name_buffer[j] <- paste(name_buffer[j-1], "/", colnames(X)[j])
      } 
    }
    DNAME <- name_buffer[length(name_buffer)]
  } else {
    DNAME<-paste("Y",  "/", "X")
  }       
  structure(list(statistic = STAT, alternative = alternative,
                 p.value = PVAL, method = METHOD, data.name = DNAME),
            class = "htest")
  
}





generate_null_model<-function(i,lr_df,X_orig,inull_model,iseed=1,irobust=FALSE,nu=5, pci_opt_method="jp"){
  set.seed(iseed+i)
  if(nrow(lr_df)>1){
    n_length = nrow(X_orig)
    null_model_parametric <-rpci(n =  n_length, beta = lr_df$beta_n, sigma_C = apply(diff(zoo::coredata(X_orig)),2,sd),rho = lr_df$rho_n[1],
                                 sigma_M =lr_df$sigma_M_n[1], sigma_R = lr_df$sigma_R_n[1],robust=irobust, nu=nu)

  } else {
    n_length = length(X_orig)
    null_model_parametric <-rpci(n =  n_length, beta = lr_df$beta_n, sigma_C = sd(diff(zoo::coredata(X_orig))),rho = lr_df$rho_n[1],
                                 sigma_M =lr_df$sigma_M_n[1], sigma_R = lr_df$sigma_R_n[1],robust=irobust, nu=nu)
    
  }
 

  X_null<-null_model_parametric[,2:ncol(null_model_parametric)]
  Y_null<-null_model_parametric[,1]
  lr_df_null<- likelihood_ratio.pci_internal(Y=Y_null,X=X_null,null_model = inull_model,pci_opt_method=pci_opt_method, robust=irobust, nu=nu)
  lr_value_null<- lr_df_null$likelihood_ratio[1]
  lr_value_null
}







