### Trying out different models using toy datasets ###
library(eforensics);

## model
## -----
model_bl    = 'bl'
model_bbl = 'bbl'
model_rn = 'rn'
model_qbl = 'qbl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=100, nCov=1, nCov.fraud=1, model=model_bl)
data = sim_data$data

## mcmc parameters
## ---------------
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2)

## samples
## -------
samples_bl = eforensics(w ~.-a, a ~ .-w,  data=data, model=model_bl, eligible.voters = "N", mcmc=mcmc);

print(summary(samples_bl)); 
summary(samples_bl, join.chains=T); 

samples_bbl = eforensics(w ~.-a, a ~ .-w,  data=data, model=model_bbl, eligible.voters = "N", mcmc=mcmc);

print(summary(samples_bbl)); 
summary(samples_bbl, join.chains=T); 

#samples_rn = eforensics(w ~.-a, a ~ .-w,  data=data, model=model_rn, eligible.voters = "N", mcmc=mcmc);

#print(summary(samples_rn)); 
#summary(samples_rn, join.chains=T);

samples_qbl = eforensics(w ~.-a, a ~ .-w,  data=data, model=model_qbl, eligible.voters = "N", mcmc=mcmc);

print(summary(samples_qbl)); 
summary(samples_qbl, join.chains=T);

