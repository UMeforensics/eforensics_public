library(eforensics);
## model
## -----
model    = 'bl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=7000, nCov=1, model=model)
data     = sim_data$data

## mcmc parameters
## ---------------
mcmc    = list(burnin=1, n.adapt=100, n.iter=100, n.chains=2)

## samples
## -------
help(eforensics)
samples    = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc)

summary(samples)
summary(samples, join.chains=T)