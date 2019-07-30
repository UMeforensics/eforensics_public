library(eforensics);
model = 'bl'

sim_data = ef_simulateData(n=7000, nCov=5, model=model)
data = sim_data$data
mcmc = list(burn.in=100, n.adapt=100, n.iter=1000, n.chains=2)
samples = eforensics(w ~.-a, a ~ .-w,  data=data, model=model, mcmc=mcmc)