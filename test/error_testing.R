### Error Testing ###
library(eforensics);

### Missing Covariates
model = 'bl'
sim_data = ef_simulateData(n=2000, nCov=1, nCov.fraud=1, model=model);
data = sim_data$data;
missing_rows <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$x1.a[missing_rows] <- NA;
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
print(try(eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc, eligible.voters="N")));


### Missing Winning Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, nCov.fraud=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$w[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
print(try(eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc, eligible.voters="N")));


### Missing Abstention Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, nCov.fraud=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$a[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
print(try(eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc, eligible.voters="N")));


### Missing Total Observations
model = 'bl' 
sim_data = ef_simulateData(n=2000, nCov=1, nCov.fraud=1, model=model); 
data = sim_data$data; 
missing_obvs <- sample(c(1:2000), size = 50, replace = FALSE, prob = NULL); 
data$N[missing_obvs] <- NA; 
mcmc = list(burn.in=1, n.adapt=100, n.iter=100, n.chains=2);
print(try(eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc, eligible.voters="N")));
