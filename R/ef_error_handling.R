
ef_check_jags <- function()
{
    if (!runjags::testjags(silent=T)$JAGS.found) {
        stop("\n\nThe election forensics package current uses JAGS, but is was not found in the machine you are running the code. Please install it and run the code again.\n\n")
    }
    invisible()
}

check_mcmc <- function(mcmc)
{
    if ("n.chains" %in% names(mcmc) & mcmc$n.chains<1) {stop("\n\nNumber of chains (n.chains) must be larger than 1. \n\n")}
    if (!"n.chains" %in% names(mcmc)) {stop("\n\nYou must include number of chains (n.chains) in the list provided in the mcmc parameter of the function eforensics(). It must must larger than one \n\n")}
    if (!"n.iter"   %in% names(mcmc)) {stop("\n\nYou must include number of iterations (n.iter) in the list provided in the mcmc parameter of the function eforensics() \n\n")}
    if (!"burn.in"  %in% names(mcmc)) {stop("\n\nYou must include the burn-in iterations (burn.in) in the list provided in the mcmc parameter of the function eforensics() \n\n")}
    if (!"n.adapt"  %in% names(mcmc)) {stop("\n\nYou must include number of adaptative steps (n.adapt) in the list provided in the mcmc parameter of the function eforensics() \n\n")}
}

check_local_fraud <- function(samples)
{
    parameters.monitored = samples.par %>% colnames(.)  %>% stringr::str_replace(string=., pattern="\\[[0-9]*]", replacement="") %>% unique
    parameters.required  = ef_get_parameters_to_monitor(model) 
    if(!all(parameters.required %in% parameters.monitored))
        stop(paste0("\n\nYou must monitor all parameters to compute the distribuion of frauds at the observation level. Use parameters='all' in the function eforensics()\n\n") )
}

