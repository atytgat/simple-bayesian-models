# Two simple Bayesian model: simulated and real data


In the first model, we compare our own implementation of a Gibbs sampler with Jags' implementation. For this purpose, we consider a linear relationship between the response and the explanatory random variables. We begin by deriving the model parameters' distributions, then, using both our own and Jags' implementation, we sample values from those distributions and compare the results with the true parameters values. Unsurprisingly, the Jags method obtains values closer to the true ones. 
The code is located in the file 'Bayesian model - simulated data.R'.

In the second, we perform a linear regression to relate the brain weight to the body weight of several animals, this time entirely with Jags' tools. The model is presented in the Report file, along with the explanations of the procedure and the analysis of the results.


