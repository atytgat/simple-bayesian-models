
##################################################################
#       Q2 e) : Gibbs algorithm and Gelman-Rubin diagnostic      #
##################################################################
rm(list=ls())

#######################
#    GIBBS SAMPLING   #
#######################

# parameters initialization
set.seed(1234)
nu = 3; n = 500
y = 170 + 5*rt(n, df=nu)
M = 1000    # number of samples


# compute the mean of the rv beta
beta.mean <- function(y.samples, lambda.samples) {
  y.samples %*% lambda.samples / sum(lambda.samples)
}

# compute the variance of the rv beta
beta.var <- function(lambda.samples, tau.sample) {
  ( tau.sample * sum(lambda.samples) )^-1
}

# compute the rate (parameter b of the gamma distribtion) of rv tau
tau.rate <- function(beta.sample, y.samples, lambda.samples) {
  ( y.samples - beta.sample )^2 %*%  lambda.samples/2 
}



#####################
#    FIRST CHAIN    #
#####################

# parameters initialization 
beta_1 = tau_1 = rep(0, M)
lbds_1 = matrix(nrow = M, ncol = n)


# initial values 
beta_1[1] = 140 
tau_1[1] = 3
lbds_1[1,] = rep(0.5, n)


# Gibbs sampling
for (i in 2:M) {
  
  # new samples lambdas
  for (j in 1:n){
    lbds_1[i,j] = rgamma( 1, (nu+1)/2, tau_1[i-1]/2 * (y[j]-beta_1[i-1])^2 + nu/2 )
  }
  
  # new sample beta
  beta_1[i] = rnorm( 1, beta.mean(y, lbds_1[i,]), beta.var(lbds_1[i,], tau_1[i-1])  )
  
  # new sample tau
  tau_1[i] = rgamma( 1, n/2, tau.rate(beta_1[i], y, lbds_1[i,]) )
  
}



#####################
#    SECOND CHAIN   #
#####################

# parameters initialization 
beta_2 = tau_2 = rep(0, M)
lbds_2 = matrix(nrow = M, ncol = n)


# initial values 
beta_2[1] = 200 
tau_2[1] = 7
lbds_2[1,] = rep(2, n)


# Gibbs sampling
for (i in 2:M) {
  
  # new samples lambdas
  for (j in 1:n){
    lbds_2[i,j] = rgamma( 1, (nu+1)/2, tau_2[i-1]/2 * (y[j]-beta_2[i-1])^2 + nu/2 )
  }
  
  # new sample beta
  beta_2[i] = rnorm( 1, beta.mean(y, lbds_2[i,]), beta.var(lbds_2[i,], tau_2[i-1])  )
  
  # new sample tau
  tau_2[i] = rgamma( 1, n/2, tau.rate(beta_2[i], y, lbds_2[i,]) )
  
}



#######################
#      TRACE PLOT     #
#######################

# trace plot for beta samples
plot(1:M, beta_1, type="l", lwd=2, col="green", xlab="iter", ylab="beta", ylim=range(140,200),  main="beta gibbs samples")
lines(1:M, beta_2, type="l", col="red", lwd=2, lty=2)
legend("topright",c("first chain","second chain"), lwd=c(2,2), col=c("green","red"), box.lty=1, lty=1:5, cex=0.8, y.intersp=0.2)

# trace plot for tau samples
plot(1:M, tau_1, type="l", lwd=2, col="green", xlab="iter", ylab="tau", ylim=range(-1,7), main="tau gibbs samples")
lines(1:M, tau_2, type="l", col="red", lwd=2, lty=2)
legend("topright",c("first chain","second chain"), lwd=c(2,2), col=c("green","red"), box.lty=1, lty=1:5, cex=0.8, y.intersp=0.2)



########################################
#       GELMAN-RUBIN DIAGNOSTIC        #
########################################
# start at iteration T < M ?


gelman.diagnostic <- function(beta) {
  # Input : parameter matrix size (I,J)
  # Output : R hat value
  
  I = dim(beta)[1]
  J = dim(beta)[2]
  
  # compute the mean of the J chains over the I iterations
  beta.mean.J = rep(0, J)
  for (j in 1:J) {
    beta.mean.J[j] = mean(beta[,j])
  }
  
  # compute mean over all indices (sequences and iterations)
  beta.mean = mean(beta.mean.J)
  
  # Between-sequence estimator
  B = I/(J-1) * sum( (beta.mean.J - beta.mean)^2 )
  
  # compute s_j^2 (variance of chain j)
  s.square.J = rep(0, J)
  for (j in 1:J) {
    s.square.J[j] = 1/(I-1) * sum( (beta[,j] - beta.mean.J[j])^2 )
  }
  
  # Within-sequence estimator
  W = 1/J * sum(s.square.J)
  
  # Unbiased estimate of the overall variance
  V.hat = (I-1)/I * W + 1/I * B
  
  # R hat
  R.hat = V.hat / W
  
}


J = 2  # number of chains
I = M  # number of iterations in each chain

beta = matrix(c(beta_1, beta_2), nrow = I) # column j = iterations of chain j
tau = matrix(c(tau_1, tau_2), nrow = I)    # column j = iterations of chain j

# Gelman-Rubin diagnostic of parameters beta and tau 
# Necessary condition for convergence : R close to 1
R.beta = gelman.diagnostic(beta)
R.beta  # 1.006316

R.tau = gelman.diagnostic(tau)
R.tau   # 0.9992421


# Gelman diagnostic plots
gelman.plot(  list(mcmc(beta_1), mcmc(beta_2) )  )
gelman.plot(  list(mcmc(tau_1), mcmc(tau_2) )  )



### Point estimation and 95% CI given respectively by the quantiles 0.5, 0.025 and 0.975 over the two chains
quantile(c(beta_1, beta_2), p = c(0.025, 0.5, 0.975))
quantile(c(tau_1, tau_2), p = c(0.025, 0.5, 0.975))
mean(c(beta_1, beta_2))
mean(c(tau_1, tau_2))


############################################################################
#       Q2 f) : Gibbs algorithm and Gelman-Rubin diagnostic using JAGS     #
############################################################################
require(R2WinBUGS)
require(coda)
require(rjags)

rm(list=ls())

# Model creation
model2f <- function(){
  for (i  in  1:n){
    y[i] ~ dnorm(beta0, lambda[i]*tau )
    lambda[i] ~ dgamma(a,b)
  }
  beta0 ~ dnorm(0,1e-6)   # approx. uniform prior over [0,1]
  tau ~ dgamma(1e-6,1e-6) # approx. uniform prior over [0,tau]
}

model2f.file <- "Model2f.bug"
write.model(model2f, model2f.file)

# Data initialization
set.seed(1234)
nu = 3; n = 500
y = 170 + 5*rt(n, df=nu)

data <- list(y = y, n = n, a = nu/2, b = nu/2)

# Starting values of the two chains
inits <- list( list(beta0 = 140, tau = 3, lambda = rep(0.5, n)), list(beta0 = 200, tau = 7, lambda = rep(2, n)) )

# JAGS model
myjags.model.2f <- jags.model(file = model2f.file, 
                              data = data, 
                              inits = inits, 
                              n.chains = length(inits), 
                              n.adapt = 1000)

# Sampling parameters beta0, tau and lambdas
samples.2f <- coda.samples(model = myjags.model.2f,
                           variable.names = c("beta0", "tau"),
                           n.iter = 10000)

# chain 1
chain1 = samples.2f[1]
beta0.chain1 = as.matrix(chain1[,'beta0'])
tau.chain1 = as.matrix(chain1[,'tau'])

# chain 2
chain2 = samples.2f[2]
beta0.chain2 = as.matrix(chain2[,'beta0'])
tau.chain2 = as.matrix(chain2[,'tau'])

# Gelman diagnostic
gelman.diag( list(mcmc(beta0.chain1), mcmc(beta0.chain2) ) )
gelman.diag( list(mcmc(tau.chain1), mcmc(tau.chain2) ) )

gelman.plot(  list(mcmc(beta0.chain1), mcmc(beta0.chain2) )  )
gelman.plot(  list(mcmc(tau.chain1), mcmc(tau.chain2) )  )


# trace plot for beta samples
plot(1:10000, beta0.chain1, type="l", lwd=2, col="green", xlab="iter", ylab="beta", ylim=range(155,190),  main="beta gibbs samples")
lines(1:10000, beta0.chain2, type="l", col="red", lwd=2, lty=2)
legend("topright",c("first chain","second chain"), lwd=c(2,2), col=c("green","red"), box.lty=1, lty=1:5, cex=0.8, y.intersp=0.2)

# trace plot for tau samples
plot(1:10000, tau.chain1, type="l", lwd=2, col="green", xlab="iter", ylab="tau", ylim=range(-0.5,2.5), main="tau gibbs samples")
lines(1:10000, tau.chain2, type="l", col="red", lwd=2, lty=2)
legend("topright",c("first chain","second chain"), lwd=c(2,2), col=c("green","red"), box.lty=1, lty=1:5, cex=0.8, y.intersp=0.2)



# geweke diagnostic (not asked)
geweke.diag(mcmc(beta0.chain1))
geweke.diag(mcmc(beta0.chain2))

geweke.diag(mcmc(tau.chain1))     
geweke.diag(mcmc(tau.chain2))

# geweke diagnostic plot 
geweke.plot(mcmc(beta0.chain1), nbins = 20)
geweke.plot(mcmc(beta0.chain2), nbins = 20)

geweke.plot(mcmc(tau.chain1), nbins = 20)
geweke.plot(mcmc(tau.chain2), nbins = 20)



### Point estimation and 95% CI given respectively by the mean and the 2.5% and 97.5% quantiles
summary(samples.2f)
