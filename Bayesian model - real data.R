
############################################################################
#                Q3 a) + b) : sampling from the posterior                  #
############################################################################
require(R2WinBUGS)
require(coda)
require(rjags)

rm(list=ls())

# loading data (/!\ change to LOCAL directory)
setwd('C:/Users/alext/OneDrive - UCL/UCL/DATA M1/Q2/LSTAT2130 - Introduction to Bayesian statistics/Projet Bayesian')
brain.data = read.table('brain.txt', header = TRUE, dec =".")
x = log(brain.data[,'bodyweight'] )   
y = log(brain.data[,'brainweight'] )


# Model creation
model3a <- function(){
  for (i  in  1:n){
    y[i] ~ dnorm(beta0 + beta1 * (x[i] - x.mean), lambda[i]*tau )
    lambda[i] ~ dgamma(a,b)
  }
  beta0 ~ dnorm(0,1e-6)   # approx. uniform prior over [0,1]
  beta1 ~ dnorm(0,1e-6)   # approx. uniform prior over [0,1]
  tau ~ dgamma(1e-6,1e-6) # approx. uniform prior over [0,tau]
}

model3a.file <- "Model3a.bug"
write.model(model3a, model3a.file)

# Data initialization
n = length(y)
nu = 3
data <- list(y = y, n = n, x = x, x.mean = mean(x), a = nu/2, b = nu/2)

# Starting values of the  chain
inits <- list( list(beta0 = 3, beta1 = 0, tau = 0.4, lambda = rep(0.5, n)) )

# JAGS model
myjags.model.3a <- jags.model(file = model3a.file, 
                              data = data, 
                              inits = inits, 
                              n.chains = length(inits), 
                              n.adapt = 1000)

# Sampling parameters beta0, tau and lambdas
samples.3a <- coda.samples(model = myjags.model.3a,
                           variable.names = c("beta0", "beta1", "tau"),
                           n.iter = 200000)
summary(samples.3a)

# MCMC chain
mcmc.chain = samples.3a[1]
beta0 = as.matrix(mcmc.chain[,'beta0'])
beta1 = as.matrix(mcmc.chain[,'beta1'])
tau = as.matrix(mcmc.chain[,'tau'])

# trace plot of the chain for each parameters
plot(1:50000, beta0, type="l", lwd=2, col="black", xlab="iter", ylab="beta0", ylim=range(3,6),  main="beta0 gibbs samples")
plot(1:50000, beta1, type="l", lwd=2, col="black", xlab="iter", ylab="beta1", ylim=range(0.2,1),  main="beta1 gibbs samples")
plot(1:50000, tau, type="l", lwd=2, col="black", xlab="iter", ylab="tau", ylim=range(-1,9), main="tau gibbs samples")



# geweke diagnostic 
geweke.diag(mcmc(beta0))
geweke.diag(mcmc(beta1))
geweke.diag(mcmc(tau))         

# geweke diagnostic plot
geweke.plot(mcmc(beta0), nbins = 20)
geweke.plot(mcmc(beta1), nbins = 20)
geweke.plot(mcmc(tau), nbins = 20)


############################################################################
#                     Q3 c) : 95% HDP computation                          #
############################################################################
# Points estimation of beta0, beta1 and tau respectively
mean(beta0) # 4.78236
mean(beta1) # 0.7089617
mean(tau)   # 1.915994

# 95% HPD
obj = mcmc.chain
HPDinterval(obj = obj, prob = 0.95)


############################################################################
#                  Q3 d) : regression lines estimation                     #
############################################################################
require(R2WinBUGS)
require(coda)
require(rjags)
require(Metrics)

rm(list=ls())

# loading data (/!\ change to LOCAL directory)
setwd('C:/Users/alext/OneDrive - UCL/UCL/DATA M1/Q2/LSTAT2130 - Introduction to Bayesian statistics/Projet Bayesian')
brain.data = read.table('brain.txt', header = TRUE, dec =".")
x = log(brain.data[,'bodyweight'] )   
y = log(brain.data[,'brainweight'] )

# Interval x for the conditional mean desired
#int = seq(min(x),max(x),by=0.5) ; K=length(int)
int = x; K=length(int)


postmean.estimate <- function(nu) {
  # input : parameter nu value
  # output : prediction of the conditional mean for all values x in int
  
  
  # Model creation
  model3d <- function(){
    for (i  in  1:n){
      y[i] ~ dnorm(beta0 + beta1 * (x[i] - x.mean), lambda[i]*tau )
      lambda[i] ~ dgamma(a,b)
    }
    
    # Credible region of interest for the conditional mean
    for (k in 1:K){
      mu.k[k] <- beta0 + beta1 * (int[k] - mean(x))
    }
    
    beta0 ~ dnorm(0,1e-6)   # approx. uniform prior over [0,1]
    beta1 ~ dnorm(0,1e-6)   # approx. uniform prior over [0,1]
    tau ~ dgamma(1e-6,1e-6) # approx. uniform prior over [0,tau]
  }
  
  model3d.file <- "Model3d.bug"
  write.model(model3d, model3d.file)
  
  # Data initialization
  n = length(y)
  nu = nu
  data <- list(y = y, n = n, x = x, x.mean = mean(x), a = nu/2, b = nu/2, int = int, K = K)
  
  # Starting values of the chain
  inits <- list( list(beta0 = 3, beta1 = 0, tau = 0.4, lambda = rep(0.5, n)) )
  
  
  # JAGS model
  myjags.model.3d <- jags.model(file = model3d.file, 
                                data = data, 
                                inits = inits, 
                                n.chains = length(inits), 
                                n.adapt = 1000)
  
  # Sampling parameters beta0, beta1, tau and mu[k]
  samples.3d <- coda.samples(model = myjags.model.3d,
                             variable.names = c("beta0", "beta1", "tau", "mu.k"),
                             n.iter = 50000)
  
  # point esimations of mu for all x values of interest 
  mcmc.chain = samples.3d[1]
  mu = rep(0, K)
  for (k in 1:K) {
    mu[k] = mean(as.matrix(mcmc.chain[,paste("mu.k[", k, "]", sep="")]))
  }
  
  return(mu)
}

### regression lines estimated with with nu = 3,5,10,15,25 in a bayesian setting
predict1 = postmean.estimate(3)
predict2 = postmean.estimate(5)
predict3 = postmean.estimate(10)
predict4 = postmean.estimate(15)
predict5 = postmean.estimate(25)

### regression line estimated in a classical frequentist setting
brain.df = log(brain.data[,3:4])
x.centered =  brain.df[,'bodyweight'] - mean(brain.df[,'bodyweight'])
brain.df[,'bodyweight'] = x.centered
lm.fit = lm(formula = brainweight~bodyweight, data = brain.df)

# prediction using the fitted model
beta0 = lm.fit$coefficients[1]
beta1 = lm.fit$coefficients[2]
predict.fr = beta0 + beta1 * (int - mean(int))

# plot of the regression lines
plot(x,y, main="linear regression lines over the brain dataset")
legend("bottomright",c("frequentist","nu=3","nu=5","nu=10","nu=15","nu=25"), lwd=c(2,2,2,2,2,2), col=c("black","blue","green","yellow","orange","red"), box.lty=1, cex=0.8, y.intersp=0.2)

lines(int, predict1, col="blue",lwd=2)
lines(int, predict2, col="green",lwd=2)
lines(int, predict3, col="yellow",lwd=2)
lines(int, predict4, col="orange",lwd=2)
lines(int, predict5, col="red",lwd=2)
lines(int, predict.fr, col="black", lwd=2)

# rmse of the regressions
rmse(y, predict.fr)
rmse(y, predict1)
rmse(y, predict2)
rmse(y, predict3)
rmse(y, predict4)
rmse(y, predict5)


# mse of the regressions
mse(y, predict.fr)
mse(y, predict1)
mse(y, predict2)
mse(y, predict3)
mse(y, predict4)
mse(y, predict5)


# mae of the regressions
mae(y, predict.fr)
mae(y, predict1)
mae(y, predict2)
mae(y, predict3)
mae(y, predict4)
mae(y, predict5)


