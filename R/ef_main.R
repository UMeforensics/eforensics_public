

## This function is used to order the formulas 3 to 6 based on the depencend variable. User can provide any order. For instance:
## eforensics(w ~1, a ~, mu.chi.s ~ 1)
## or
## eforensics(w ~1, a ~, mu.iota.m ~ 1)
## mu.chi.s must be formula6, and mu.iota.m must be formula3, but in these two examples, both become formula3 because the user is not specifying formula3=mu.iota.m~1 in the second example nor formula6=mu.chi.s ~ 1 in the second
## this function use the dependent variables (mu.chi.s, mu.chi.m, mu.iota.s, mu.iota.m) to order the formulas
## It returns a list of four formulas ordered correctly
order.formulas <- function(formula3, formula4, formula5, formula6)
{
    if (is.null(formula3)) {formula3=NA}
    if (is.null(formula4)) {formula4=NA}
    if (is.null(formula5)) {formula5=NA}
    if (is.null(formula6)) {formula6=NA}
    formulas = list(formula3, formula4, formula5, formula6)
    idx = purrr::map_dbl(.x=formulas, function(.x)  dplyr::case_when(stringr::str_detect(.x, pattern="mu.iota.m") ~ 1,
                                                                     stringr::str_detect(.x, pattern="mu.iota.s") ~ 2,
                                                                     stringr::str_detect(.x, pattern="mu.chi.m")  ~ 3,
                                                                     stringr::str_detect(.x, pattern="mu.chi.s")  ~ 4,
                                                                     all(is.na(.x %>% as.character)) ~ 999 ## formula NULL
                                                                     ) %>%
                                                    .[!is.na(.)]
                         )
    formulas.final = list(NA,NA,NA,NA)
    for (i in 1:length(idx))
    {
        if (idx[i]!=999) {
            formulas.final[[idx[i]]] = formulas[[i]]
        }
    }
    if ( is.na(formulas.final[[1]]) ) {formulas.final[[1]] = mu.iota.m ~ 1}
    if ( is.na(formulas.final[[2]]) ) {formulas.final[[2]] = mu.iota.s ~ 1}
    if ( is.na(formulas.final[[3]]) ) {formulas.final[[3]] = mu.chi.m ~ 1}
    if ( is.na(formulas.final[[4]]) ) {formulas.final[[4]] = mu.chi.s ~ 1}

    return(formulas.final)
}

## convergence diagnostics
## -----------------------
ef_get_diagnostic <- function(samples, mcmc.conv.diagnostic, mcmc.conv.parameters, mcmcse.conv.precision, mcmcse.combine)
{
    options(warn=-1)
    on.exit(options(warn=0))
    if (mcmc.conv.diagnostic == "PSRF") {
        summ              = summary(samples) %>% base::data.frame(Parameter=rownames(.), .) %>% tibble::as_tibble()
        results           = summ %>% dplyr::select(Parameter, psrf, Mean)
        converged         = all(results$psrf < 1.05)
    }
    if (mcmc.conv.diagnostic == "MCMCSE") {
        conv.parameters.regexp = paste0(mcmc.conv.parameters, collapse="|")
        ## combine the chains and compute MCMC SE
        ## It passes the diagnostic if each parameter in the combined chains passes the criterium
        if (mcmcse.combine) {
            samples                = runjags::combine.mcmc(samples)
            ## mcmcse                 = mcmcse::mcse.multi(samples)$est ## posterior average
            mcmcse                 = mcmcse::mcse.multi(samples)$cov %>% diag(.) %>% sqrt(.)
            results                = mcmcse[stringr::str_subset(names(mcmcse), pattern=conv.parameters.regexp)] %>%
                tibble::enframe(., name=c("Parameter")) %>%
                dplyr::rename(MCMCSE=value) %>%
                dplyr::mutate(MCMCSE.criterium = mcmcse.conv.precision,
                              Converged        = MCMCSE < mcmcse.conv.precision)
        }else{
            ## compute MCMC for each chain, get the maximum MCMC SE for each parameter and chain, and compare the maximum with the criterium
            ## It passes the diagnostic if each parameter in each chain passes the criterium
            parameters = samples$mcmc[[1]] %>% colnames
            ncols      = length(samples$mcmc)
            results = samples$mcmc %>%
                purrr::map(.x=., function(.x, .y) mcmcse::mcse.multi(.x)$cov %>% diag(.) %>% sqrt(.) %>% tibble::as_tibble(.) ) %>%
                dplyr::bind_cols(.) %>%
                dplyr::rename_all(list(~paste0("chain.", 1:ncols) )) %>%
                dplyr::mutate(MCMCSE.max       = pmax(!!!rlang::syms(names(.)), na.rm=TRUE),
                              MCMCSE.criterium = mcmcse.conv.precision,
                              Converged        = MCMCSE.max < mcmcse.conv.precision,
                              Parameter = parameters) %>%
                dplyr::select(Parameter, dplyr::everything())  %>%
                dplyr::filter(stringr::str_detect(Parameter, pattern=conv.parameters.regexp))
        }
        converged = all(results$Converged)
    }
    return(list(diagnostic=mcmc.conv.diagnostic, results = results, converged=converged))
}

ef_print_diagnostic <- function(diagnostic)
{
    msg     = paste0('\n','Convergence diagnostic: ', diagnostic$diagnostic, '\n'); cat(msg)
    print(diagnostic$results)
}
get_Z <- function(samples)
{
    ## replace matrix Z (sample size) x (number of interation) to a single column matrix with z.hat
    ## the estimated cluster of Zi. Zi is classified in the cluster it has highest estimated
    ## posterior probability to belong to
    for (i in 1:length(samples))
    {
        z            = samples[[i]][,base::grepl(pattern='Z.[0-9]*.', x=colnames(samples[[i]]))]
        k.hat        =   base::apply(z, 2, function(zi) which.max(c(sum(zi==1), sum(zi==2), sum(zi==3)) ) )
        piZi         = t(base::apply(z, 2, function(zi) c(sum(zi==1), sum(zi==2), sum(zi==3))/length(zi)))
        colnames(piZi )= c("pi[Zi1]", "pi[Zi2]", "pi[Zi3]")
        samp         = samples[[i]][,!base::grepl(pattern='Z.[0-9]*.', x=colnames(samples[[i]]))]
        samp         = list(parameters=coda::as.mcmc(samp), k.hat=k.hat, piZi=piZi)
        samples[[i]] = samp
    }
    return(samples)
}

get_Z_qbl <- function(samplez)
{
  ## replace matrix Z (sample size) x (number of interation) to a single column matrix with z.hat
  ## the estimated cluster of Zi. Zi is classified in the cluster it has highest estimated
  ## posterior probability to belong to
  #Create holders for individual proportions of fraud
  samplez2 <- list()
  im <- matrix(ncol = dim(samplez[[1]][,base::grepl(pattern='^iota[.]m.[0-9]*.', x=colnames(samplez[[1]]))])[2], nrow = 0)
  is <- matrix(ncol = dim(samplez[[1]][,base::grepl(pattern='^iota[.]s.[0-9]*.', x=colnames(samplez[[1]]))])[2], nrow = 0)
  cm <- matrix(ncol = dim(samplez[[1]][,base::grepl(pattern='^chi[.]m.[0-9]*.', x=colnames(samplez[[1]]))])[2], nrow = 0)
  cs <- matrix(ncol = dim(samplez[[1]][,base::grepl(pattern='^chi[.]s.[0-9]*.', x=colnames(samplez[[1]]))])[2], nrow = 0)
  zz <- matrix(ncol = dim(samplez[[1]][,base::grepl(pattern='Z.[0-9]*.', x=colnames(samplez[[1]]))])[2], nrow = 0)
  for (i in 1:length(samplez))
  {
    z            = samplez[[i]][,base::grepl(pattern='Z.[0-9]*.', x=colnames(samplez[[i]]))]
    im           = rbind(im,samplez[[i]][,base::grepl(pattern='^iota[.]m.[0-9]*.', x=colnames(samplez[[i]]))])
    is           = rbind(is,samplez[[i]][,base::grepl(pattern='^iota[.]s.[0-9]*.', x=colnames(samplez[[i]]))])
    cm           = rbind(cm,samplez[[i]][,base::grepl(pattern='^chi[.]m.[0-9]*.', x=colnames(samplez[[i]]))])
    cs           = rbind(cs,samplez[[i]][,base::grepl(pattern='^chi[.]s.[0-9]*.', x=colnames(samplez[[i]]))])
    zz           = rbind(zz,samplez[[i]][,base::grepl(pattern='Z.[0-9]*.', x=colnames(samplez[[i]]))])
    k.hat        =   base::apply(z, 2, function(zi) which.max(c(sum(zi==1), sum(zi==2), sum(zi==3)) ) )
    piZi         = t(base::apply(z, 2, function(zi) c(sum(zi==1), sum(zi==2), sum(zi==3))/length(zi)))
    colnames(piZi )= c("pi[Zi1]", "pi[Zi2]", "pi[Zi3]")
    samp         = samplez[[i]][,!base::grepl(pattern='Z.[0-9]*.|^iota[.]m.[0-9]*.|^iota[.]s.[0-9]*.|^chi[.]m.[0-9]*.|^chi[.]s.[0-9]*.', x=colnames(samplez[[i]]))]
    samp         = list(parameters=coda::as.mcmc(samp), k.hat=k.hat, piZi=piZi)
    samplez2[[i]] = samp
  }
  #Get appropriate fraud props
  prop.fraud.m <- zz
  prop.fraud.s <- zz
  for(i in 1:dim(zz)[1]){
    for(j in 1:dim(zz)[2]){
      zzi <- zz[i,j]
      if(zzi == 1){
        prop.fraud.m[i,j] <- 0
        prop.fraud.s[i,j] <- 0
      }else{
        if(zzi == 2){
          prop.fraud.m[i,j] <- im[i,j]
          prop.fraud.s[i,j] <- is[i,j]
        }else{
          prop.fraud.m[i,j] <- cm[i,j]
          prop.fraud.s[i,j] <- cs[i,j]
        }
      }
    }
  }
  get.hpd.top <- function(x){
    return(coda::HPDinterval(coda::as.mcmc(x), .95)[2])
  }
  get.hpd.low <- function(x){
    return(coda::HPDinterval(coda::as.mcmc(x), .95)[1])
  }
  #Summarize the output
  #Posterior mean, 95%HPD, and quantiles .01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99
  post.mean.m <- base::apply(prop.fraud.m,2,mean)
  post.hpdu.m <- base::apply(prop.fraud.m,2,get.hpd.top)
  post.hpdl.m <- base::apply(prop.fraud.m,2,get.hpd.low)
  post.quan.m <- base::apply(prop.fraud.m,2,stats::quantile,c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99))
  post.mean.s <- base::apply(prop.fraud.s,2,mean)
  post.hpdu.s <- base::apply(prop.fraud.s,2,get.hpd.top)
  post.hpdl.s <- base::apply(prop.fraud.s,2,get.hpd.low)
  post.quan.s <- base::apply(prop.fraud.s,2,stats::quantile,c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99))
  #Put everything together for presentation
  post.m <- cbind(post.mean.m,post.hpdl.m,post.hpdu.m,t(post.quan.m))
  rownames(post.m) <- as.character(seq(1,dim(post.m)[1]))
  colnames(post.m) <- c("Mean","95%HPDLow","95%HPDHigh","Q1%","Q2.5%","Q5%","Q10%","Q25%","Q50%","Q75%","Q90%","Q95%","Q97.5%","Q99%")
  post.s <- cbind(post.mean.s,post.hpdl.s,post.hpdu.s,t(post.quan.s))
  rownames(post.s) <- as.character(seq(1,dim(post.s)[1]))
  colnames(post.s) <- c("Mean","95%HPDLow","95%HPDHigh","Q1%","Q2.5%","Q5%","Q10%","Q25%","Q50%","Q75%","Q90%","Q95%","Q97.5%","Q99%")
  frauds <- list()
  frauds$Manufactured <- post.m
  frauds$Stolen <- post.s
  out.list <- list()
  out.list$samples <- samplez2
  out.list$frauds <- frauds
  return(out.list)
}

create_list <- function(samples)
{
    for (i in 1:length(samples))
    {
        samples[[i]] = list(parameters=samples[[i]])
    }
    return(samples)
}

ef_get_parameters_to_monitor <- function(model, all=FALSE)
{
    ## Restricted normal models
    ## ------------------------
    ## if(model == 'rn')                  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn')                  parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s", "alpha")
    if(model == 'rn_no_alpha')         parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
    ## if(model == 'rn_no_alpha')         parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s")
    if(model == 'normal')              parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")

    if(model == 'rn_sep')              parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'normal_sep')          parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_alpha_sep')     parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s")

    if(model == 'rn_no_scaled')        parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    if(model == 'rn_no_scaled_sep')    parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")

    ## binomial models
    ## ---------------
    ## if(model == 'bl')                  parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s", "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")
    if(model == 'bl')               parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
    if(model == 'bl.rd')               parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")

    ## overdispersion model (beta binomial)
    ## --------------------
    ## if(model == 'bbl')              parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s","nf.var")
    if(model == 'bbl')              parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
    if(model == 'bbl.rd')           parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")
    if(model == 'bl')               parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
    if(model == 'qbl')              parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")

    ## varying dimension models
    ## ------------------------
    if(model == 'bl.vd')               parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s","psi.i")

    if(model == 'rn.vd')               parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_scaled.vd')     parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'normal.vd')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_alpha.vd')      parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s", "psi.i")

    if(model == 'rn_sep.vd')           parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_scaled_sep.vd') parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'normal_sep.vd')       parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha", "psi.i")
    if(model == 'rn_no_alpha_sep.vd')  parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "mu.iota.s", "mu.chi.s", "sigma.iota.s", "psi.i")

    if(all) parameters = c(parameters, 'Z')

    return(parameters)
}

ef_get_parameters_to_monitor_c <- function(model)
{
    if(model == 'rn'){
        parameters = c('pi', 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "sigma.iota.m", "sigma.tau", "sigma.nu", "alpha")
    }else{
        if(model == 'bl'){
            parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
        }else{
            if(model == "bbl"){
                parameters = c("pi", 'beta.tau', 'beta.nu', "beta.iota.m", "beta.iota.s", "beta.chi.m", "beta.chi.s")
                ## parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")
            }else{
                if(model == "bl.rd"){
                    parameters = c("pi", 'beta.tau', 'beta.nu', "mu.iota.m",  "mu.chi.m", "mu.iota.s", "mu.chi.s")
                }else{
                    parameters = c("pi")
                }
            }
        }
    }
    return(parameters)
}

get_model <- function(model)
{
    if (model == 'rn')           return(rn())
    if (model == 'rn_no_scaled') return(rn_no_scaled())
    if (model == 'normal')       return(normal())
    if (model == 'rn_no_alpha')  return(rn_no_alpha())

    if (model == 'rn_sep')           return(rn_sep())
    if (model == 'rn_no_scaled_sep') return(rn_no_scaled_sep())
    if (model == 'normal_sep')       return(normal_sep())
    if (model == 'rn_no_alpha_sep')  return(rn_no_alpha_sep())

    if (model == 'bl')          return(bl())
    if (model == 'bl.rd')          return(bl.rd())
    if (model == 'bl_fc')       return(bl_cov())

    if (model == 'bbl')          return(bbl())
    if (model == 'qbl')          return(qbl())
    if (model == 'bbl.rd')       return(bbl.rd())

    if (model == 'bl.vd')       return(bl.vd())

    if (model == 'rn.vd')           return(rn.vd())
    if (model == 'rn_no_scaled.vd') return(rn_no_scaled.vd())
    if (model == 'normal.vd')       return(normal.vd())
    if (model == 'rn_no_alpha.vd')  return(rn_no_alpha.vd())

    if (model == 'rn_sep.vd')           return(rn_sep.vd())
    if (model == 'rn_no_scaled_sep.vd') return(rn_no_scaled_sep.vd())
    if (model == 'normal_sep.vd')       return(normal_sep.vd())
    if (model == 'rn_no_alpha_sep.vd')  return(rn_no_alpha_sep.vd())

}

getRegMatrix <- function(func.call, data, weights, formula_number=1)
{
    args <- names(func.call)
    ## creating the dependent variable and the covariates matrix from the fomula 1
    f = paste0('formula', formula_number, sep='')
    idx.args  <- match(c(f,  "data", "weights"), args, 0L)
    func.call <- func.call[c(1L, idx.args)]
    names(func.call)[names(func.call)==f] = "formula"
    func.call$drop.unused.levels <- TRUE
    func.call[[1L]] <- quote(stats::model.frame)
    func.call[[3]] = quote(data)
    reg.matrix <- eval(func.call, parent.frame())
    ## response variable
    y   <- stats::model.response(reg.matrix, "numeric")
    ## weights
    w   <- as.vector(stats::model.weights(reg.matrix))
    if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
    offset <- as.vector(stats::model.offset(func.call))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y))
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", length(offset), NROW(y)), domain = NA)
    }
    ## covariates
    mt1    <- attr(reg.matrix, "terms")
    if (stats::is.empty.model(mt1)) {
        x <- matrix(1, ncol=1,nrow=nrow(y))
        results <- list(coefficients = if (is.matrix(y)) matrix(, 0, 3) else numeric(), residuals = y, fitted.values = 0 * y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            results$fitted.values <- offset
            results$residuals <- y - offset
        }
    } else {
        x <- stats::model.matrix(mt1, reg.matrix, contrasts)
    }
    return(list(y=y, X=x, w=w))
}


ef_parameters_to_check_convergence <- function(mcmc.conv.parameters)
{
    currently_implemented = c("pi")
    ## currently_implemented = c("pi", "alpha")
	if (!all(mcmc.conv.parameters %in% currently_implemented)) {
        warning(paste0("\n\n\n *** Note about automatic convergence checks *** \n Current implementation check convergence only for the following parameters: ", paste0(currently_implemented, collapse=", "), "\n\n"  ) )
    }
}

## {{{ docs }}}

#' Election Forensics Finite Mixture Model
#'
#' This function estimates a finite mixture model of election fraud
#'
#'
#' @param formula1 an object of the class \code{formula} as used in \code{\link{lm}}. The dependent variable of this formula must the number (counts) or proportion of votes for the party or candidate that won the election. If counts are used, the model must be from the binomial family (see \code{model} parameter below). If proportions are provided, the model must be from the normal family (see \code{model} parameter below)
#' @param formula2 an object of the class \code{formula} as used in \code{\link{lm}}. The dependent variable of this formula must the number (counts) or proportion of abstention.  The type (count or proportion) must be the same as the independent variable in \code{formula1}
#' @param formula3 See description below
#' @param formula4 See description below
#' @param formula5 See description below
#' @param formula6 See description below
#' \describe{
#'   \item{Formulas 3 to 6}{There are four other possible formulas to use: formula3, formula4, formula5, formula6}
#'   \item{formula3}{an object of the class \code{formula} as used in \code{\link{lm}}. The left-hand side (LHS) of the formula must be mu.iota.m (see example). The mu.iota.m is the probability of incremental fraud by manufacturing votes and it is a latent variable in the model. By specifying the LHS with that variable, the functional automatically identifies that formula as formula3. Default is \code{NULL} and it means that probability is not affected by election unit (ballot box, polling place, etc) covariate}
#'   \item{formula4}{an object of the class \code{formula} as used in \code{\link{lm}}. The left-hand side (LHS) of the formula must be mu.iota.s (see example). The mu.iota.s is the probability of incremental fraud by stealing votes from the opposition and it is a latent variable in the model. By specifying the LHS with that variable, the functional automatically identifies that formula as formula4. Default is \code{NULL} and it means that probability is not affected by election unit (ballot box, polling place, etc) covariate}
#'   \item{formula5}{an object of the class \code{formula} as used in \code{\link{lm}}. The left-hand side (LHS) of the formula must be mu.chi.m (see example). The mu.chi.m is the probability of extreme fraud by manufacturing votes and it is a latent variable in the model. By specifying the LHS with that variable, the functional automatically identifies that formula as formula5. Default is \code{NULL} and it means that probability is not affected by election unit (ballot box, polling place, etc) covariate}
#'   \item{formula6}{an object of the class \code{formula} as used in \code{\link{lm}}. The left-hand side (LHS) of the formula must be mu.chi.s (see example). The mu.chi.s is the probability of extreme fraud by stealing votes from the opposition and it is a latent variable in the model. By specifying the LHS with that variable, the functional automatically identifies that formula as formula6. Default is \code{NULL} and it means that probability is not affected by election unit (ballot box, polling place, etc) covariate}
#' }
#' @param model a string with the model ID to use in the estimation. There are three current choices: \code{qbl}, \code{bl}, and \code{rn}.  \code{qbl} is the default and recommended choice.  For a description of each model, see \code{ef_models_desc}. 
#' @param data a data.frame with the independent variables (voters for the winner and abstention) and the covariates. If the independent variables are counts, the it is necessary to provide the total number of eligible voters (see parameter \code{eligible.voters})
#' @param eligible.voters string with the name of the variable in the data that contains the number of eligible voters. Default is \code{NULL}, but it is required if the independent variables (voters for the winner and abstention) are counts
#' @param weights Deprecated.
#' @param mcmc a list containing \code{n.iter}, which is the number of iterations for the MCMC, \code{burn.in} for the burn-in period of the MCMC chain, \code{n.adapt} indicating the number of adaptative steps before the estimation (see \code{rjags}), and \code{n.chains}, an integer indicating the number of chains to use (default 1).
#' @param parameters a string vector with the names of the parameters to monitor. When \code{NULL}, it will monitor all the parameters, except the Z's. When \code{parameters='all'} (default), it will monitor all parameters, including Z, which is necessary to classify the observations as fraudulent cases or not.
#' @param na.action Deprecated.
#' @param get.dic Deprecated.
#' @param parComp Logical.  If \code{parComp = TRUE}, then chains are computed in parallel using the runjags parallel method.  This opens \code{n.chains} instances of JAGS.  In practice, a max of 4 unique chains can be run due to the way in which JAGS generates initial values.  If \code{parComp = FALSE}, chains are run sequentially using the runjags interruptible method.
#' @param autoConv Logical.  If \code{autoConv = TRUE}, chains are run until convergence criteria are met.  Currently, chains are run for a single period equal to \code{burn.in} iterations and monitored for \code{n.iter} iterations.  If \code{mcmc.conv.diagnostic = "MCMCSE"}, MCMCSE values are calculated for each parameter in \code{mcmc.conv.diagnostic}.  If all values are less than \code{mcmcse.conv.precision} then the chain is stopped and the chain is run for \code{n.iter} more iterations monitoring all values specified by \code{parameters}.  If the MCMCSE for any parameter is higher than \code{mcmcse.conv.precision}, then the chain is run for \code{burn.in} + \code{n.iter} more iterations and the MCMCSE is again checked.  This is repeated, at most, \code{max.auto} times.  If the MCMCSE condition is not met by \code{max.auto} attempts, a warning message is printed and the chains are run \code{n.iter} more times with all parameters monitored.  If \code{mcmc.conv.diagnostic = "PSRF"}, the same procedure occurs checking that all PSRF values are less than 1.05.
#' @param max.auto Integer.  Number of subsequent tries to achieve the convergence conditions outlined by \code{autoConv}.  After \code{max.auto} failures, a warning is thrown and the chain is run \code{n.iter} more times monitoring all specified parameters.
#' @param mcmc.conv.diagnostic a string with the method to use to evaluate convergence. Currenctly, \code{PSRF} and \code{MCMCSE} (default) are implemented.
#' @param mcmc.conv.parameters string vector with the name of the parameters to check for convergence using the MCMC standard error. Default is \code{pi},
#' @param mcmcse.conv.precision numeric, the value of the precision criterion to evaluate convergence using the MCMC standard error. The MCMC std. error of all parameters included in \code{mcmcse.conv.parameter} must be below the threshold defined by the value of \code{mcmc.conv.precision} (default is 0.05) to pass the convergence diagnostic.
#' @param mcmcse.combine boolean, if \code{TRUE}, the MCMCSE is computed after the chains are combined. Otherwise, the MCMC std. error is computed for each chain, and the maximum std. error of each parameter is used for the diagnostic
#'
#' @return The function returns a nested list of class \code{eforensics} with length equal to the number of chains. Each sublist contains up to three named objects:
#' \describe{
#'   \item{parameters}{A \code{mcmc} object that contains the posterior draws for all monitored parameters except for the individual fraud classifications.}
#'   \item{k.hat}{A vector that contains the posterior modal classification for each observation.  1 corresponds to no fraud, 2 corresponds to incremental fraud, and 3 corresponds to extreme fraud.}
#'   \item{piZi}{A matrix with three columns that contains the posterior probability of belonging to each class for each observation.}
#' }
#'   
#' If \code{model = "qbl"} or \code{model = "bl"}, the proportion of frauds estimated at each observation is returned.  These values can be accessed for object \code{foo} using \code{attr(foo,"frauds")}.  This attribute is a two element list that contains the estimated proportion of votes that are Stolen and Manufactured.  Posterior means, HPD intervals, and posterior quantiles are returned for each observation in the data set.  These quantities are automatically aggregated over all chains. 
#'
#' @references
#' Flegal, J. M., Haran, M., & Jones, G. L., Markov chain monte carlo: can we trust the third significant figure?, Statistical Science, 23(2), 250–260 (2008).
#' Brooks, S. P., & Gelman, A., General methods for monitoring convergence of iterative simulations, Journal of computational and graphical statistics, 7(4), 434–455 (1998).
#' @examples
#' 
#' set.seed(12345)
#' library(eforensics)
#' model    = 'qbl'
#' 
#' ## simulate data
#' ## -------------
#' set.seed(12345)
#' sim_data = ef_simulateData(n=250, nCov=1, nCov.fraud=1,
#'               model="bbl", overdispersion = 100, pi = c(.95,.04,.01))
#' data     = sim_data$data
#' 
#' ## mcmc parameters
#' ## ---------------
#' mcmc    = list(burn.in=1000, n.adapt=1000, n.iter=1000, n.chains=2)
#' 
#' ## samples
#' ## -------
#' ## help(eforensics)
#' 
#' samples    = eforensics(
#'   w ~ x1.w ,
#'   a ~ x1.a,
#'   mu.iota.m ~ x1.iota.m,
#'   mu.iota.s ~ x1.iota.s,
#'   mu.chi.m  ~ x1.chi.m,
#'   mu.chi.s  ~ x1.chi.s,
#'   data=data,
#'   eligible.voters="N",
#'   model="qbl",
#'   mcmc=mcmc,
#'   parameters = "all",
#'   parComp = TRUE,
#'   autoConv = TRUE,
#'   max.auto = 10,
#'   mcmc.conv.diagnostic = "MCMCSE",
#'   mcmc.conv.parameters = c("pi"),
#'   mcmcse.conv.precision = .05,
#'   mcmcse.combine = FALSE
#' )
#' 
#' #Summaries for each of the monitored parameters
#' #Look at each chain separately
#' summary(samples)
#' #Combine the chains
#' summary(samples, join.chains=T)
#' 
#' #Look at the estimated fraud proportions for each observation
#' attr(samples,"frauds")
#' #Look at Manufactured and Stolen separately
#' attr(samples,"frauds")$Manufactured
#' attr(samples,"frauds")$Stolen
#' 
#' #How accurate is the classification?
#' 
#' #Get the true categories
#' true_z <- sim_data$latent$z
#' 
#' #What is the modal estimate for the class?
#' num_z <- (samples[[1]]$piZi*1000) + (samples[[2]]$piZi*1000)
#' max_z <- apply(num_z,1,which.max)
#' 
#' #How accurate is the modal classification?
#' table(true_z,max_z)
#' 
#' #How accurately do we uncover the proportion of frauds for each observation?
#' 
#' #Manufactured
#' true_man <- ((true_z == 1)*0) + ((true_z == 2)*sim_data$latent$iota.m) + 
#'             ((true_z == 3)*sim_data$latent$chi.m)
#' 
#' #What is the posterior mean proportion of manufactured votes
#' pred_man <- attr(samples,"frauds")$Manufactured[,1]
#' 
#' #Are they close?
#' plot(true_man, pred_man, xlab = "True Proportion Manufactured Votes", 
#'      ylab = "Estimated Proportion Manufactured Votes")
#' 
#' #Stolen
#' true_stolen <- ((true_z == 1)*0) + ((true_z == 2)*sim_data$latent$iota.s) + 
#'                ((true_z == 3)*sim_data$latent$chi.s)
#' 
#' #What is the posterior mean proportion of manufactured votes
#' pred_stolen <- attr(samples,"frauds")$Stolen[,1]
#' 
#' #Are they close?
#' plot(true_stolen, pred_stolen, xlab = "True Proportion Stolen Votes", 
#'      ylab = "Estimated Proportion Stolen Votes")
#'
#' @export

## }}}
eforensics   <- function(formula1, formula2, formula3=NULL, formula4=NULL, formula5=NULL, formula6=NULL, data, eligible.voters=NULL, weights=NULL, mcmc, model = "qbl",
                         parameters="all", na.action="exclude", get.dic = 1000, parComp = TRUE, autoConv = TRUE, max.auto = 10,
                         mcmc.conv.diagnostic="MCMCSE", mcmc.conv.parameters = "pi", mcmcse.conv.precision = 0.05, mcmcse.combine=FALSE)
{

    ## error handling
    check_mcmc(mcmc)
    check_nchains_for_psrf(mcmc, mcmc.conv.diagnostic)

    ## check parameters used to check convergence
    ef_parameters_to_check_convergence(mcmc.conv.parameters)
    
    options(warn=-1)
    on.exit(options(warn=0))
    ## check if JAGS is installed
    ef_check_jags()

    ## order the formulas (see comment in the order.formulas() function)
    formulas = order.formulas(formula3, formula4, formula5, formula6)
    formula3 = formulas[[1]] %>% stats::as.formula(.)
    formula4 = formulas[[2]] %>% stats::as.formula(.)
    formula5 = formulas[[3]] %>% stats::as.formula(.)
    formula6 = formulas[[4]] %>% stats::as.formula(.)

    ## create a placeholder for weights it if is not provided
    if (is.null(weights)) {data = data %>% dplyr::mutate(weights = 1)}

    #if(parComp == TRUE){
        eforensics_main_par(formula1, formula2, formula3, formula4, formula5, formula6, data, eligible.voters, weights, mcmc, model, parameters, na.action, get.dic, autoConv, max.auto,
                            mcmc.conv.diagnostic=mcmc.conv.diagnostic, mcmc.conv.parameters=mcmc.conv.parameters, mcmcse.conv.precision=mcmcse.conv.precision, mcmcse.combine=mcmcse.combine, parComp = parComp)
    # }else{
    #     eforensics_main(formula1, formula2, formula3, formula4, formula5, formula6, data, eligible.voters, weights, mcmc, model, parameters, na.action, get.dic,
    #                     mcmc.conv.diagnostic=mcmc.conv.diagnostic, mcmc.conv.parameters=mcmc.conv.parameters, mcmcse.conv.precision=mcmcse.conv.precision, mcmcse.combine=mcmcse.combine)
    # }

}

#Dear future Kevin, this is a branch of the code that used to use rjags instead of rjags to implement the eforensics functions.  The other branch is superior and can do parallel computation with little backend work.

# eforensics_main   <- function(formula1, formula2, formula3, formula4, formula5, formula6, data, eligible.voters=NULL, weights, mcmc, model, parameters=NULL, na.action="exclude", get.dic = 1000,
#                               mcmc.conv.diagnostic, mcmc.conv.parameters, mcmcse.conv.precision, mcmcse.combine=mcmcse.combine)
# {
#     ## error handling
#     check_mcmc(mcmc)
#     options(warn=-1)
#     on.exit(options(warn=0))
#     ## check if JAGS is installed
#     ef_check_jags()
# 
#     ## constructing formulas 3 to 6 placeholder in the data for the latent variables mu.iota.s, mu.iota.m, mu.chi.s, mu.chi.m. This is needed to construct the design matrix
#     data$mu.iota.m = 1
#     data$mu.iota.s = 1
#     data$mu.chi.m = 1
#     data$mu.chi.s = 1
# 
#     ## ## construct the regression matrices (data.frames) based on the formula provided
#     ## ## -----------------------------------------------------------------------------
#     func.call <- match.call(expand.dots = FALSE)
#     ## votes for the winner
#     mat     = getRegMatrix(func.call, data, weights, formula_number=1)
#     w       = mat$y
#     Xw      = mat$X
#     weightw = mat$w
#     ## abstention
#     mat     = getRegMatrix(func.call, data, weights, formula_number=2)
#     a       = mat$y
#     Xa      = mat$X
#     weighta = mat$w
# 
#     ## incremental fraud manufactures (mu.iota.m)
#     mat      = getRegMatrix(func.call, data, weights, formula_number=3)
#     X.iota.m = mat$X
#     weighta  = mat$w
#     ## incremental fraud stolen (mu.iota.s)
#     mat      = getRegMatrix(func.call, data, weights, formula_number=4)
#     X.iota.s = mat$X
#     weighta  = mat$w
#     ## incremental fraud manufactures (mu.chi.m)
#     mat      = getRegMatrix(func.call, data, weights, formula_number=5)
#     X.chi.m = mat$X
#     weighta  = mat$w
#     ## incremental fraud stolen (mu.chi.s)
#     mat      = getRegMatrix(func.call, data, weights, formula_number=6)
#     X.chi.s = mat$X
#     weighta  = mat$w
#     dat    = list(w = w, a = a,
#                   Xa       = as.matrix(Xa)      , dxa       = ncol(Xa),
#                   Xw       = as.matrix(Xw)      , dxw       = ncol(Xw),
#                   X.iota.m = as.matrix(X.iota.m), dx.iota.m = ncol(X.iota.m),
#                   X.iota.s = as.matrix(X.iota.s), dx.iota.s = ncol(X.iota.s),
#                   X.chi.m  = as.matrix(X.chi.m),  dx.chi.m  = ncol(X.chi.m),
#                   X.chi.s  = as.matrix(X.chi.s),  dx.chi.s  = ncol(X.chi.s),
#                   n = length(w))
#     if(!is.null(eligible.voters)){
#         data = data %>% dplyr::rename(eligible.voters = !!eligible.voters)
#         dat$N = data$eligible.voters
#     }else{
#         ## check if model use counts, and require elebigle voters
#         ## ------------------------------------------------------
#         if (stringr::str_detect(model, pattern="bl")) {
#             stop("\nThe parameter 'eligible.voters' must be provided for models based on binomial distributions\n\n")
#         }
#     }
#     data = dat
#     ## get parameters to monitor
#     ## -------------------------
#     if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)
#     if(parameters[1] == 'all') parameters = ef_get_parameters_to_monitor(model, all=TRUE)
# 
#     ## get model
#     ## ---------
#     model.name = model
#     model      = get_model(model.name)
# 
#     ## Debug/Monitoring message --------------------------
#     msg <- paste0('\n','Burn-in: ', mcmc$burn.in, '\n'); cat(msg)
#     ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
#     msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
#     msg <- paste0('\n','MCMC in progress ....', '\n'); cat(msg)
#     ## ---------------------------------------------------
# 
#     ## MCMC
#     ## ----
#     time.init    = Sys.time()
#     cat('\nCompiling the model...\n')     ; sim = rjags::jags.model(file=textConnection(model), data = data, n.adapt=mcmc$n.adapt, n.chain=mcmc$n.chains);
#     cat('\nUpdating MCMC (burn-in) ...\n'); stats::update(sim, n.iter = mcmc$burn.in)
#     cat('\nDrawing the samples...\n')     ; samples = rjags::coda.samples(model=sim, variable.names=parameters, n.iter=mcmc$n.iter)
#     T.mcmc = Sys.time() - time.init
#     if(get.dic != 0){
#         cat('\nDrawing DIC samples...\n') ; dic.samples = rjags::dic.samples(model=sim, n.iter=get.dic, type = "popt")
#     }else{
#         dic.samples = NULL
#     }
#     T.mcmc = Sys.time() - time.init
# 
#     if(!is.null(parameters) & "Z" %in% parameters)
#         samples = get_Z(samples)
#     else
#         samples = create_list(samples)
#     class(samples) = "eforensics"
# 
#     attr(samples, "formula.w") = formula1
#     attr(samples, "formula.a") = formula2
#     attr(samples, "model")     = model.name
#     if (model.name %in% c("rn")) {
#         attr(samples, "terms")     = c("alpha", colnames(X.chi.m), colnames(X.iota.m), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
#     }else{
#         attr(samples, "terms")     = c(colnames(X.chi.m), colnames(X.chi.s), colnames(X.iota.m), colnames(X.iota.s), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
#     }
#     attr(samples, "dic")       = dic.samples
#     ## include alpha in terms name if using rn model
# 
#     cat("\n\nEstimation Completed\n\n")
#     return(samples)
# }

eforensics_main_par   <- function(formula1, formula2, formula3, formula4, formula5, formula6, data, eligible.voters=NULL, weights, mcmc, model, parameters=NULL, na.action="exclude", get.dic = 1000, autoConv = TRUE, max.auto = 10,
                                  mcmc.conv.diagnostic, mcmc.conv.parameters, mcmcse.conv.precision, mcmcse.combine=mcmcse.combine, parComp)
{
    ## error handling
    check_mcmc(mcmc)
    options(warn=-1)
    on.exit(options(warn=0))
    ## check if JAGS is installed
    ef_check_jags()

    ## constructing formulas 3 to 6 placeholder in the data for the latent variables mu.iota.s, mu.iota.m, mu.chi.s, mu.chi.m. This is needed to construct the design matrix
    data$mu.iota.m = 1
    data$mu.iota.s = 1
    data$mu.chi.m = 1
    data$mu.chi.s = 1

    ## ## construct the regression matrices (data.frames) based on the formula provided
    ## ## -----------------------------------------------------------------------------
    func.call <- match.call(expand.dots = FALSE)
    ## votes for the winner
    mat     = getRegMatrix(func.call, data, weights, formula_number=1)
    w       = mat$y
    Xw      = mat$X
    weightw = mat$w
    ## abstention
    mat     = getRegMatrix(func.call, data, weights, formula_number=2)
    a       = mat$y
    Xa      = mat$X
    weighta = mat$w

    ## incremental fraud manufactures (mu.iota.m)
    mat      = getRegMatrix(func.call, data, weights, formula_number=3)
    X.iota.m = mat$X
    weighta  = mat$w
    ## incremental fraud stolen (mu.iota.s)
    mat      = getRegMatrix(func.call, data, weights, formula_number=4)
    X.iota.s = mat$X
    weighta  = mat$w
    ## incremental fraud manufactures (mu.chi.m)
    mat      = getRegMatrix(func.call, data, weights, formula_number=5)
    X.chi.m = mat$X
    weighta  = mat$w
    ## incremental fraud stolen (mu.chi.s)
    mat      = getRegMatrix(func.call, data, weights, formula_number=6)
    X.chi.s = mat$X
    weighta  = mat$w
    dat    = list(w = w, a = a,
                  Xa       = as.matrix(Xa)      , dxa       = ncol(Xa),
                  Xw       = as.matrix(Xw)      , dxw       = ncol(Xw),
                  X.iota.m = as.matrix(X.iota.m), dx.iota.m = ncol(X.iota.m),
                  X.iota.s = as.matrix(X.iota.s), dx.iota.s = ncol(X.iota.s),
                  X.chi.m  = as.matrix(X.chi.m),  dx.chi.m  = ncol(X.chi.m),
                  X.chi.s  = as.matrix(X.chi.s),  dx.chi.s  = ncol(X.chi.s),
                  n = length(w))
    if(!is.null(eligible.voters)){
        data = data %>% dplyr::rename(eligible.voters = !!eligible.voters)
        dat$N = data$eligible.voters
    }else{
        ## check if model use counts, and require elebigle voters
        ## ------------------------------------------------------
        if (stringr::str_detect(model, pattern="bl")) {
            stop("\nThe parameter 'eligible.voters' must be provided for models based on binomial distributions\n\n")
        }
    }
    data = dat
    ## get parameters to monitor
    ## -------------------------
    if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)
    if(parameters[1] == 'all') parameters = ef_get_parameters_to_monitor(model, all=TRUE)
    ## parameters_c = ef_get_parameters_to_monitor_c(model)

    ## get model
    ## ---------
    model.name = model
    model      = get_model(model.name)

    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n','Burn-in: ', mcmc$burn.in, '\n'); cat(msg)
    ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
    msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
    msg <- paste0('\n','MCMC in progress ....', '\n\n'); cat(msg)
    ## ---------------------------------------------------
    
    ##Use parComp flag to determine runjags method
    if(parComp == TRUE){
      rjMethod <- "parallel"
    }else{
      rjMethod <- "interruptible"
    }

    ## MCMC
    ## ----
    time.init    = Sys.time()
    ## Run for a prespecified amount of time
    runjags::runjags.options(inits.warning=FALSE, rng.warning=FALSE)
    ## manke sure max.auto is not negative of a fraction. If fraction, get the smaller integer larger than the fraction
    max.auto = max(0, ceiling(max.auto))
    if(autoConv == TRUE){
        trial      = 0
        presim     = runjags::run.jags(model = model, monitor = mcmc.conv.parameters, data = data, n.chains = mcmc$n.chains, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = .5, method = rjMethod)
        diagnostic = ef_get_diagnostic(presim, mcmc.conv.diagnostic, mcmc.conv.parameters, mcmcse.conv.precision, mcmcse.combine=mcmcse.combine)
        while(trial < max.auto & !diagnostic$converged){
            ## Messages about convergence
            ef_print_diagnostic(diagnostic)
            trial      = trial + 1
            cat(paste("\nConvergence diagnostic requirement not met. Extending the chains (attempt ", trial, " of ", max.auto, ") ... \n\n",sep = ""))
            ## extending the chain
            presim = runjags::extend.jags(presim, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = 10, method = rjMethod)
            diagnostic = ef_get_diagnostic(presim, mcmc.conv.diagnostic, mcmc.conv.parameters, mcmcse.conv.precision, mcmcse.combine=mcmcse.combine)
        }
        ## Final sample to keep
        cat("\nBurnin Finished.\nCapturing the samples ...\n\n")
        if(model.name %in% c("qbl","bl")){
          #Add individual model params
          amp <- c(parameters,"iota.m","iota.s","chi.m","chi.s")
        }else{
          amp <- parameters
        }
        samples = runjags::extend.jags(presim, add.monitor = c(amp), burnin = 0, sample = mcmc$n.iter, adapt = mcmc$n.adapt, thin = 1, method = rjMethod, jags.refresh = 10)
        ## Print final diagnostic
        ef_print_diagnostic(diagnostic)
        if(!diagnostic$converged){
            cat("\n\n")
            cat(paste("NOTE:\nConvergence diagnostics indicate that the chain(s) didn\'t converge.\n", sep=""))
            cat(paste0("Use these estimated results with caution and think about re-running eforensics with more chains, samples, or a longer burnin.", sep = ""))
            cat("\n")
        }
    }else{
      ## sample to keep
      if(model.name %in% c("qbl","bl")){
        #Add individual model params
        amp <- c(parameters,"iota.m","iota.s","chi.m","chi.s")
      }else{
        amp <- parameters
      }
      samples = runjags::run.jags(model = model, monitor = amp, data = data, n.chains = mcmc$n.chains, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = .5, method = rjMethod)
      cat("\n\n")
      cat(paste("NOTE: Convergence diagnostics were not computed'. Use results with caution.\n", sep=""))
      cat(paste0("Set 'autoConv=T' to compute diagnostic automatically. See help(eforensics).", sep = ""))
      cat("\n")
    }
    dic.samples = NULL
    T.mcmc = Sys.time() - time.init

    if(!is.null(parameters) & "Z" %in% parameters & model.name != "qbl" & model.name != "bl"){
      samples = get_Z(samples[[1]])
    }else{
      if(!is.null(parameters) & "Z" %in% parameters & model.name %in% c("qbl","bl")){
        samples = get_Z_qbl(samplez = samples[[1]])
        frauds = samples$frauds
        samples = samples$samples
      }else{
        samples = create_list(samples[[1]])
      }
    }
    class(samples) = "eforensics"
    
    
    attr(samples, "formula.w") = formula1
    attr(samples, "formula.a") = formula2
    attr(samples, "model")     = model.name
    if (model.name %in% c("rn")) {
      attr(samples, "terms")     = c("alpha", colnames(X.chi.m), colnames(X.iota.m), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
    }else{
      ## attr(samples, "terms")     = c(colnames(X.chi.m), colnames(X.chi.s), colnames(X.iota.m), colnames(X.iota.s), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
      attr(samples, "terms")     = c("No Fraud", "Incremental Fraud", "Extreme Fraud",
                                     colnames(Xw), colnames(Xa),
                                     colnames(X.iota.m), colnames(X.chi.m),
                                     colnames(X.chi.m),  colnames(X.chi.s)    )
    }
    if(model.name %in% c("qbl","bl")){
      attr(samples,"frauds") = frauds
    }else{
      flist <- list()
      flist$Manufactured <- NULL
      flist$Stolen <- NULL
      attr(samples,"frauds") = flist
    }

    ## attr(samples, "formula.w") = formula1
    ## attr(samples, "formula.a") = formula2
    ## attr(samples, "model")     = model.name
    ## attr(samples, "terms")     = c(colnames(X.chi.m), colnames(X.chi.s), colnames(X.iota.m), colnames(X.iota.s), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
    ## attr(samples, "dic")       = "No DIC can be calculated using parallel chains"

    cat("\n\nEstimation Completed\n\n")
    return(samples)
}


## {{{ parallel (old) }}}

## eforensics_par   <- function(formula1, formula2, data, weights, mcmc, model, parameters=NULL, na.action="exclude")
## {

##     ## check if JAGS is installed
##     ef_check_jags()

##     ## ## construct the regression matrices (data.frames) based on the formula provided
##     ## ## -----------------------------------------------------------------------------
##     func.call <- match.call(expand.dots = FALSE)
##     mat     = getRegMatrix(func.call, data, weights, formula_number=1)
##     w       = mat$y
##     Xw      = mat$X
##     weightw = mat$w
##     mat     = getRegMatrix(func.call, data, weights, formula_number=2)
##     a       = mat$y
##     Xa      = mat$X
##     weighta = mat$w
##     if(model == 'bl'){
##         data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w), N = data$N)
##     }else{
##         data    = list(w = w, a = a, Xa = as.matrix(Xa), Xw = as.matrix(Xw), dxw = ncol(Xw), dxa = ncol(Xa), n = length(w))
##     }


##     ## get parameters to monitor
##     ## -------------------------
##     if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)

##     ## get model
##     ## ---------
##     model = get_model(model)

##     ## Debug/Monitoring message --------------------------
##     msg <- paste0('\n','Burn-in: ', mcmc$burn.in, '\n'); cat(msg)
##     ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
##     msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
##     msg <- paste0('\n','MCMC in progress ....', '\n'); cat(msg)
##     ## ---------------------------------------------------

##     ## MCMC (parallel)
##     ## ----
##     cl <- parallel::makePSOCKcluster(min(parallel::detectCores()-1, mcmc$n.chains))
##     samples = dclone::jags.parfit(cl       = cl,
##                                   data     = data,
##                                   params   = parameters,
##                                   model    = textConnection(model),
##                                   n.adapt  = mcmc$n.adapt,
##                                   n.chains = mcmc$n.chains,
##                                   n.update = mcmc$burn.in,
##                                   n.iter   = mcmc$n.iter,
##                                   thin     = 1)
##     parallel::stopCluster(cl)


##     ## computing summary
##     ## -----------------
##     ## summary <- list(summary = summary(samples), HPD = coda::HPDinterval(samples))
##     ## results = list(samples=samples, stat=summary, time.elapsed=T.mcmc)
##     results = samples

##     cat("\n\nEstimation Completed\n\n")
##     return(results)
## }

## }}}
