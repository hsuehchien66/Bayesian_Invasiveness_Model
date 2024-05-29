library(rstan)
library(shinystan)
setwd("/Users/hc14/Documents/PhD_project/Invasiveness/Stan_Bayesian/simulation")

# generating fake data
N1 <- 100
Y1 <- rnorm(N, 1.6, 0.2) # Y is height
hist(Y1)

# load model
model <- stan_model("simulation.stan")

# pass the data to stan model
options(mc.cores=4) # assign multiple cores
fit <- sampling(model, list(N=N1, Y=Y1), iter=200, chains=4)

# inspect the results: posterior distribution of the parameters
params <- extract(fit)
hist(params$mu)

# use shinystan to get more details
launch_shinystan(fit)



