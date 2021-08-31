# clean up workspace environment
rm(list = ls())

# load library 
library(tidyverse)
library(quantmod)
library(LambertW) # for KS test of t-distribution (reference: https://www.rdocumentation.org/packages/LambertW/versions/0.6.5/topics/ks.test.t)
library(MASS) # for cov.trob() profile likelihood func
library(rugarch) 
library(copula)
library(mnormt)
library(fGarch) # for pstd

# load data
getSymbols("^GSPC", src = "yahoo",from = as.Date("2006-01-01"), to = as.Date("2021-04-30")) # S&P500
getSymbols("^STOXX", src = "yahoo",from = as.Date("2006-01-01"), to = as.Date("2021-04-30")) # EURO STOXX
getSymbols("^HSI", src = "yahoo",from = as.Date("2006-01-01"), to = as.Date("2021-04-30")) # HongKong Hang Seng Index
getSymbols("^N225", src = "yahoo",from = as.Date("2006-01-01"), to = as.Date("2021-04-30")) # Japan Nikkei 225

total = merge.xts(GSPC$GSPC.Close, STOXX$STOXX.Close, HSI$HSI.Close, N225$N225.Close) # adjusted closing price for each index
total = na.omit(total) # omit na
return_net = na.omit(ROC(total, type = 'discrete')) # net returns
return_log = na.omit(ROC(total, type = 'continuous')) # log returns

# fit univariate distributions to each returns
# check normality by Shapiro-Wilks test
shapiro.test(as.vector(return_net$GSPC.Close)) 
shapiro.test(as.vector(return_net$STOXX.Close))
shapiro.test(as.vector(return_net$HSI.Close))
shapiro.test(as.vector(return_net$N225.Close))
# result: all reject null hypothesis (i.e. samples from normal distribution) and not follow normal distributions

# QQ plot for normal quantiles
par(mfrow = c(2,2))
name= c("S&P 500", "Euro STOXX", "Hang Seng index", "Nikkei 225")
for (i in 1:4){
  qqnorm(return_net[,i], main = paste(name[i], "QQ plot"))
  qqline(return_net[,i], col = "red", lwd = 3)
  
}

# check goodness of fit for t-distribution to each returns
ks.test.t(as.vector(return_net$GSPC.Close))
ks.test.t(as.vector(return_net$STOXX.Close))
ks.test.t(as.vector(return_net$HSI.Close))
ks.test.t(as.vector(return_net$N225.Close))
# result: all have relatively large p-value (greater than 0.05) and pass the KS test for t-distribution

# fit the std t-distribution to each returns and extract the three parameters
GSPC_t = as.numeric(fitdist(distribution = "std", return_net$GSPC.Close)$pars)
STOXX_t = as.numeric(fitdist(distribution = "std", return_net$STOXX.Close)$pars)
HSI_t = as.numeric(fitdist(distribution = "std", return_net$HSI.Close)$pars)
N225_t = as.numeric(fitdist(distribution = "std", return_net$N225.Close)$pars)

# QQ plot for t-quantiles
# S&P 500 
n = length(return_net$GSPC.Close) # number of sample points
grid = (1:n)/(n+1) 

df = GSPC_t[3] # fitted tail index v for S&P 500
qqplot(as.numeric(return_net$GSPC.Close), qt(grid, df = df), main = "S&P 500 QQ plot", xlab = 'net returns', ylab = 't-quantiles')
lm_fit = lm(qt(c(0.25,0.75), df = df)~ quantile(return_net$GSPC.Close, c(0.25,0.75))) # fit a linear model for the different quantiles between 
abline(lm_fit, col = 'red')

# Euro STOXX
n = length(return_net$STOXX.Close) # number of sample points
grid = (1:n)/(n+1) 

df = STOXX_t[3] # fitted tail index v for Euro STOXX
qqplot(as.numeric(return_net$STOXX.Close), qt(grid, df = df), main = "Euro STOXX QQ plot", xlab = 'net returns', ylab = 't-quantiles')
lm_fit = lm(qt(c(0.25,0.75), df = df)~ quantile(return_net$STOXX.Close, c(0.25,0.75))) # fit a linear model for the different quantiles between 
abline(lm_fit, col = 'red')

# HongKong Hang Seng Index
n = length(return_net$HSI.Close) # number of sample points
grid = (1:n)/(n+1) 

df = HSI_t[3] # fitted tail index v for HongKong Hang Seng
qqplot(as.numeric(return_net$HSI.Close), qt(grid, df = df), main = "Hang Seng Index QQ plot", xlab = 'net returns', ylab = 't-quantiles')
lm_fit = lm(qt(c(0.25,0.75), df = df)~ quantile(return_net$HSI.Close, c(0.25,0.75))) # fit a linear model for the different quantiles between 
abline(lm_fit, col = 'red')

# Nikkei 225
n = length(return_net$N225.Close) # number of sample points
grid = (1:n)/(n+1) 

df = N225_t[3] # fitted tail index v for Japan Nikkei 225
qqplot(as.numeric(return_net$N225.Close), qt(grid, df = df), main = "Nikkei 225 QQ plot", xlab = 'net returns', ylab = 't-quantiles')
lm_fit = lm(qt(c(0.25,0.75), df = df)~ quantile(return_net$N225.Close, c(0.25,0.75))) # fit a linear model for the different quantiles between 
abline(lm_fit, col = 'red')

par(mfrow = c(1,1))
# fit multivariate t-distribution by profile likelihood method 
df = seq(2.1, 5, 0.01) # potential degree of freedom v
n = length(df) # number of d.o.f
loglik = rep(NA,n) # the profile likelihood for each d.o.f
for (i in 1:n){
  est = cov.trob(return_net, nu = df[i])
  loglik[i] = sum(log(dmt(return_net, mean = est$center, S = est$cov, df = df[i])))
}
nuhat = df[which.max(loglik)] # the MLE of degree of freedom v
nuhat # optimal tail index and refit the profile likelihood function
est = cov.trob(return_net, nu = nuhat) # MLES of the mean vector mu and the scale matrix in fitted multivariate-t distribution
est

# plot for profile log-likelihood
plot(df,2*loglik-2*max(loglik), type = "l", xlab = "tail index v", ylab = "2*log-likelihood - max", main = "Profile log-likelihood w.r.t tail index v")
abline(v=df[which.max(loglik)], col = "red", lwd = 3)

# fit t-copula and meta-t distribution
# Apply Kendall's tau to determine the correlation matrix
Kendall_tau = cor(return_net, method = "kendall")
correlation_matrix = sin((pi/2)*Kendall_tau)
correlation_matrix_upper = c(correlation_matrix[1,2], correlation_matrix[1,3], correlation_matrix[1,4],correlation_matrix[2,3],correlation_matrix[2,4],correlation_matrix[3,4]) # store useful info from the transformed correlation matrix

# fit the t-copula
t_copula = tCopula(correlation_matrix_upper, dim = 4, dispstr = "un", df = 4) # fit t copula with an unstructured correlation matrix 
cdf_returns = cbind(pdist(distribution = "std", return_net$GSPC.Close, GSPC_t[1], GSPC_t[2], GSPC[3]),
                 pdist(distribution = "std", return_net$STOXX.Close, STOXX_t[1], STOXX_t[2], STOXX_t[3]),
                 pdist(distribution = "std", return_net$HSI.Close, HSI_t[1], HSI_t[2], HSI_t[3]),
                 pdist(distribution = "std", return_net$N225.Close, N225_t[1], N225_t[2], N225_t[3])) # transform the net returns by their estimated CDF to obtain approximately uniformly-distributed variables

fitCopula = fitCopula(t_copula, cdf_returns, method = "ml", start = c(correlation_matrix_upper,4))
fitCopula # method = 'ml' means parameters are fitted by MLE method 

# fit a meta-t distribution with an unstructured correlation matrix and std univariate marginal distributions to net returns
meta_t_dist = mvdc(t_copula, c("std","std","std","std"), list(
  list(mean = GSPC_t[1], sd = GSPC_t[2], nu = GSPC_t[3]),
  list(mean = STOXX_t[1], sd = STOXX_t[2], nu = STOXX_t[3]),
  list(mean = HSI_t[1], sd = HSI_t[2], nu = HSI_t[3]),
  list(mean = N225_t[1], sd = N225_t[2], nu = N225_t[3])
)) # defines the multivariate distribution from the t-copula and the "std" univariate marginal distributions

start = c(GSPC_t, STOXX_t, HSI_t, N225_t, fitCopula@estimate) # start with 4 fitted univariate t-distribution parameters and t-copula parameters
fit_meta_t = fitMvdc(as.matrix(return_net), meta_t_dist, start) # took a long time to execute...
signif(fit_meta_t@estimate,3) # fitted parameters of the marginal distributions and the t-copula in the order as our input
signif(sqrt(diag(fit_meta_t@var.est)),3) # standard errors of each estimated parameter
