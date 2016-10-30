
#' Asymmetric Laplace Working Likelihood For 
#' Linear Quantile Regression 
#' 
#' @param pars regression coefficient vector
#' @param y the response vector
#' @param x the design matrix with one in the first column corresponding to the intercept
#' @param tau the quantile level
#' @param sig scale parameter sigma
#' @return the working log (asymmetric Laplace) likelihood function (the part involving the regression coefficients)
#' @details The asymmetric Laplace working likelihood
#'  is proportional to exponential
#' of the negative quantile objective function for 
#' linear quantile regression
#' @export
#' @references Yang, Y., Wang, H. and He, X. (2015). Posterior inference in Bayesian quantile regression with asymmetric Laplace
#'   likelihood. International Statistical Review, 2015. doi: 10.1111/insr.12114. 

li_reg=function(pars, y, x, tau, sig){
  r = y - x%*%pars
  rhot=(tau-1*{r<0})*r
  out = -sum(rhot)/sig
  return(out)
}

#' Estimation of the Scale Parameter for Quantile Regression
#' @param y the response vector
#' @param x the design matrix with one in the first column corresponding to the intercept
#' @return the estimated scale parameter sigma 
#' @export
#' @importFrom quantreg rq
sigma_est = function(y, x)
{
  #estimate sigma
  bhat=rq(y~x-1,tau=0.5)$coef
  sig=-li_reg(bhat, y, x, tau=0.5, 1)/length(y)
  return(sig)
}

#' Asymmetric-Laplace-type Working Likelihood For Cenaored Quantile Regression
#' 
#' Asymmetric-Laplace-type working likelihood for linear 
#' quantile regression with responses 
#' subject to left censoring at zero
#' 
#' @details The asymmetric-Laplace-type working likelihood is proportional to exponential
#' of the negative Powell objective function for censored quantile regression
#' @references Powell, J. L. (1986). Censored regression quantiles. Journal of Econometrics, 32, 143-155.
#' @references Yang, Y., Wang, H. and He, X. (2015). Posterior inference in Bayesian quantile regression with asymmetric Laplace
#'   likelihood. International Statistical Review, 2015. doi: 10.1111/insr.12114. 
#' @param pars regression coefficient vector
#' @param y the response vector
#' @param x the design matrix with one in the first column corresponding to the intercept
#' @param tau the quantile level
#' @param sig scale parameter sigma
#' @return the working log (asymmetric Laplace-type) likelihood function (the part involving the regression coefficients)
#' @export
#' 
li_powell=function(pars,y,x,tau,sig){
  xb=x%*%pars
  r=y-pmax(xb,0)
  rhot=(tau-1*{r<0})*r
  return(-sum(rhot)/sig)
}
#' Estimation of the Scale Parameter for Censored Quantile Regression
#' @param y the response vector
#' @param x the design matrix with one in the first column corresponding to the intercept
#' @return the estimated scale parameter sigma 
#' @export
#' @importFrom quantreg rq
sigma_est_powell = function(y, x)
{
  n=length(y)
  bhat=crq(Curv(y,rep(0,n))~x-1,taus=0.5,method="Pow")$coef
  sig=-li_powell(bhat, y, x, tau=0.5, 1)/n
  return(sig)  
}


#' Adjusted Bayesian Quantile Regression 
#' 
#' Bayesian quantile regression 
#' based on asymmetric Laplace likelihood 
#' with posterior variance adjustment
#' @details 
#' The function returns the unadjusted and adjusted posterior standard deviation, and unadjusted and 
#' adjusted credible intervals for Bayesian quantile regression based on asymmetric Laplace working likelihood.
#' @param y the response vector
#' @param x the design matrix. If the first column of x is not all ones, a column of ones will be added.
#' @param tau the quantile level of interest
#' @param niter integer: number of iterations to run the chain for. Default 20000.
#' @param burn_in  integer: discard the first burn_in values. Default 100.
#' @param prop_cov  covariance matrix giving the covariance of the proposal distribution. 
#' This matrix need not be positive definite. 
#' If the covariance structure of the target distribution is known (approximately), it can be given here. 
#' If not given, the diagonal will be estimated via the Fisher information matrix.
#' @param level nominal confidence level for the credible interval
#' @export
#' @importFrom quantreg rq
#' @importFrom MHadaptive Metro_Hastings
#' @importFrom stats cov qnorm
#' @importFrom coda effectiveSize
#' @return A list of the following commponents is returned
#' @return estpar: posterior mean of the regression coefficient vector
#' @return PSD: posterior standard deviation without adjustment
#' @return PSD.adj: posterior standard deviation with adjustment
#' @return CI.BAL: credible interval without adjustment
#' @return CI.BAL.adj: credible interval with adjustment
#' @return sig: estimated scale parameter
#' @return MCMCsize: effective size of the chain
#' @references Yang, Y., Wang, H. and He, X. (2015). Posterior inference in Bayesian quantile regression with asymmetric Laplace
#'   likelihood. International Statistical Review, 2015. doi: 10.1111/insr.12114. 
#' @examples 
#' #A simulation example
#' library(AdjBQR)
#' n=200
#' set.seed(12368819)
#' x1 = rnorm(n)  
#' x2 = rnorm(n)  
#' y=2*x1+2*x2+rt(n,df=3)
#' x = cbind(1, x1, x2)
#' ## Bayesian quantile regression based on asymmetric Laplace likelihood
#' BQR(y, x, tau=0.5, level=0.9)
BQR = function(y, x, tau, niter=20000, burn_in=4000, prop_cov=NULL, level=0.9)
{
  x = as.matrix(x)
  
  #some input check
  stopifnot(length(y) == nrow(x))
  stopifnot(niter > burn_in)
  
  if(sum(abs(x[,1]-1))>0) x=cbind(1, x)
  p = ncol(x)
  n = length(y)
  inipars=rq(y~x-1,tau)$coef
  sig = sigma_est(y, x)
  par.names = paste("b", (1:p)-1, sep="")
  
  #run mcmc
  if(is.null(prop_cov))
    mcmc_r=Metro_Hastings(li_func=li_reg,pars=inipars,par_names=par.names,
                          iterations=niter,burn_in=burn_in,quiet=T,y=y, 
                          x=x,tau=tau,sig=sig)
  else
    mcmc_r=Metro_Hastings(li_func=li_reg,pars=inipars,
                          prop_sigma=prop_cov,
                          par_names=par.names,
                          iterations=niter,burn_in=burn_in,quiet=T,y=y, 
                          x=x,tau=tau,sig=sig)
  
  postchain=mcmc_r$trace
  
  #posterior covariance matrix (uncorrected)
  Sigma = cov(postchain)
  
  #adjusted covariance matrix
  D0 = t(x)%*%x/n
  Sigma.adj = n*tau*(1-tau)/sig^2*Sigma%*%D0%*%Sigma
  
  estpar = colMeans(postchain)
  
  #BAL Interval without adjustment
  PSD = sqrt(diag(Sigma))
  LB = estpar - qnorm(1-(1-level)/2)*PSD
  UB = estpar + qnorm(1-(1-level)/2)*PSD
  CI.BAL = cbind(LB, UB)
  
  #BAL Interval with adjustment
  PSD.adj =sqrt(diag(Sigma.adj))
  LB = estpar - qnorm(1-(1-level)/2)*PSD.adj
  UB = estpar + qnorm(1-(1-level)/2)*PSD.adj
  CI.BAL.adj = cbind(LB, UB)
  
  MCMCsize=effectiveSize(postchain)
  out = list(estpar=estpar, PSD=PSD, PSD.adj=PSD.adj, CI.BAL=CI.BAL, 
             CI.BAL.adj=CI.BAL.adj, sig=sig, MCMCsize=MCMCsize)
  return(out)
}


#' Adjusted Bayesian Censored Quantile Regression 
#' 
#' Bayesian quantile regression based on asymmetric-Laplace-type 
#' likelihood with posterior variance adjustment
#' @details 
#' The function returns the unadjusted and adjusted posterior standard deviation, and unadjusted and 
#' adjusted credible intervals for Bayesian censored quantile regression based on asymmetric-Laplace-type
#'  working likelihood. The asymmetric-Laplace-type likelihood is based on the objective function
#'  of the Powell's estimator in Powell (1986).
#' @param y the observed response vector that is left censored at zero
#' @param x the design matrix. If the first column of x is not all ones, a column of ones will be added.
#' @param tau the quantile level of interest
#' @param niter integer: number of iterations to run the chain for. Default 20000.
#' @param burn_in  integer: discard the first burn_in values. Default 100.
#' @param prop_cov  covariance matrix giving the covariance of the proposal distribution. 
#' This matrix need not be positive definite. 
#' If the covariance structure of the target distribution is known (approximately), it can be given here. 
#' If not given, the diagonal will be estimated via the Fisher information matrix.
#' @param level nominal confidence level for the credible interval
#' @export
#' @importFrom quantreg crq Curv
#' @importFrom MHadaptive Metro_Hastings
#' @importFrom stats cov qnorm
#' @importFrom coda effectiveSize
#' @importFrom survival survfit
#' @return A list of the following commponents is returned
#' @return estpar: posterior mean of the regression coefficient vector
#' @return PSD: posterior standard deviation without adjustment
#' @return PSD.adj: posterior standard deviation with adjustment
#' @return CI.BAL: credible interval without adjustment
#' @return CI.BAL.adj: credible interval with adjustment
#' @return sig: estimated scale parameter
#' @return MCMCsize: effective size of the chain
#' @references Powell, J. L. (1986). Censored regression quantiles. Journal of Econometrics, 32, 143-155.
#' @references Yang, Y., Wang, H. and He, X. (2015). Posterior inference in Bayesian quantile regression with asymmetric Laplace
#'   likelihood. International Statistical Review, 2015. doi: 10.1111/insr.12114. 
#' @examples 
#' #A simulation example
#' library(AdjBQR)
#' n=200
#' set.seed(12368819)
#' x1=rnorm(n)
#' x2=rnorm(n)
#' ystar=3/4+2*x1+3*x2+rt(n,df=3)
#' y=ystar*(ystar>0)
#' delta=1*(ystar>0)
#' x = cbind(x1, x2)
#' ## Bayesian censored quantile regression based on asymmetric-Laplace-type likelihood
#' BCQR(y, x, tau=0.5, level=0.9)

BCQR = function(y, x, tau, niter=20000, burn_in=4000, prop_cov=NULL, level=0.9)
{
  x = as.matrix(x)
  if(sum(abs(x[,1]-1))>0) x=cbind(1, x)
  p = ncol(x)
  n = length(y)
  
  #some input check
  stopifnot(length(y) == nrow(x))
  stopifnot(niter > burn_in)
  
  
  inipars=crq(Curv(y,rep(0,length(y)))~x-1,taus=tau,method="Pow")$coef
  sig = sigma_est_powell(y, x)
  par.names = paste("b", (1:p)-1, sep="")
  
  #run mcmc
  if(is.null(prop_cov))
    mcmc_r=Metro_Hastings(li_func=li_powell, pars=inipars, 
                          par_names=par.names,iterations=niter,
                          burn_in=burn_in,quiet=T,y=y, x=x, 
                          tau=tau, sig=sig) 
  else
    mcmc_r=Metro_Hastings(li_func=li_powell, pars=inipars, 
                          prop_sigma=prop_cov,
                          par_names=par.names,iterations=niter,
                          burn_in=burn_in,quiet=T,y=y, x=x, 
                          tau=tau, sig=sig)     
  postchain=mcmc_r$trace  
  estpar=colMeans(postchain)
  
  #posterior covariance matrix (uncorrected)
  Sigma = cov(postchain)
  
  #adjusted covariance matrix
  fit = x%*%estpar
  idx = which(fit>0)
  D0 = t(x[idx,])%*%x[idx,]/n
  Sigma.adj = n*tau*(1-tau)/sig^2*Sigma%*%D0%*%Sigma
  
  #BAL Interval without adjustment
  PSD = sqrt(diag(Sigma))
  LB = estpar - qnorm(1-(1-level)/2)*PSD
  UB = estpar + qnorm(1-(1-level)/2)*PSD
  CI.BAL = cbind(LB, UB)
  
  #BAL Interval with adjustment
  PSD.adj = sqrt(diag(Sigma.adj))
  LB = estpar - qnorm(1-(1-level)/2)*PSD.adj
  UB = estpar + qnorm(1-(1-level)/2)*PSD.adj
  CI.BAL.adj = cbind(LB, UB)
  
  MCMCsize=effectiveSize(postchain)
  out = list(estpar=estpar, PSD=PSD, PSD.adj=PSD.adj, CI.BAL=CI.BAL, 
             CI.BAL.adj=CI.BAL.adj, sig=sig, MCMCsize=MCMCsize)
  return(out)
}




