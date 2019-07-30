Overview
========

Instalation
===========

**NOTE:** [JAGS](https://sourceforge.net/projects/mcmc-jags/) is
required to use the current version of this package

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
# Install the development version (requires the package "devtools", so install it first if it is not installed already)
devtools::install_github("DiogoFerrari/eforensics")
# If you don't want to update the dependencies, use: (you may need to install some dependencies manually)
devtools::install_github("DiogoFerrari/eforensics", dependencies=F)
```

Usage
=====

In R, check `help(eforensics)`

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
library(eforensics)


ef_models()
model  = 'bl'

## simulating data and parameters
## ------------------------------
help(ef_simulateData)
sim_data = ef_simulateData(n=1000, nCov=1, nCov.fraud=2,  model=model)
data     = sim_data$data %>% tibble::as_data_frame() 

## mcmc parameters
## ---------------
mcmc    = list(burn.in=1, n.adapt=10, n.iter=200, n.chains=2)

## samples
## -------
help(eforensics)
samples    = eforensics(
    w ~ x1.w ,
    a ~ x1.a,
    mu.iota.s ~ x1.iota.s + x2.iota.s,
    mu.chi.m  ~ x1.chi.m + x2.chi.m ,
    mu.chi.s  ~ x1.chi.s + x2.chi.s ,
    mu.iota.m ~ x1.iota.m + x2.iota.m ,
    data=data,
    elegible.voters="N",
    model=model, mcmc=mcmc, get.dic=0,
    parameters = "all")

## summary
## -------
summary(samples)
summary(samples, join.chains=T)

## compare with the true
## ---------------------
True = ef_get_true(sim_data) 
summary(samples, join.chains=T) %>% dplyr::left_join(., True , by="Parameter") 

## plots
## -----
ef_plot(samples)
ef_plot(samples, True)
ef_plot(samples, True, parse=T)
ef_plot(samples, plots = c("Abstention and Vote",  "Fraud Probability" ))
ef_plot(samples, plots = c( "Fraud Probability" ))
```


