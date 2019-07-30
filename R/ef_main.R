

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
ef_get_diagnostic <- function(samples, diagnostic='PSRF')
{
    if (diagnostic == "PSRF") {
        summ              = summary(samples) %>% base::data.frame(Parameter=rownames(.), .) %>% tibble::as_tibble() 
        results           = summ %>% dplyr::select(Parameter, psrf, Mean) 
        converged         = all(results$psrf < 1.05)
    }
    return(list(diagnostic=diagnostic, results = results, converged=converged))
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
#' @param model a string with the model ID to use in the estimation. Run \code{\link{ef_models()}} to see the list and descriptions of the models available.
#' @param data a dara.frame with the independent variables (voters for the winner and abstention) and the covariates. If the independent variables are counts, the it is necessary to provide the total number of elegible voters (see parameter \code{elegible.voters})
#' @param elegible.voters string with the name of the variable in the data that contains the number of elegible voters. Default is \code{NULL}, but it is required if the independent variables (voters for the winner and abstention) are counts
#' @param weights (not used)
#' @param mcmc a list containing \code{n.iter}, which is the number of iterations for the MCMC, \code{burn.in} for the burn-in period of the MCMC chain, \code{n.adapt} indicating the number of adaptative steps before the estimation (see \code{\link{rjags}})
#' @param parameters a string vector with the names of the parameters to monitor. When \code{NULL}, it will monitor all the parameters, except the Z's. When \code{parameters='all'} (default), it will monitor all parameters, including Z, which is necessary to classify the observations as fraudulent cases or not.
#' @param na.action (not used)
#' @param get.dic logical.  If get.dic is FALSE, no DIC is calculated.  If get.dic is an integer greater than 0, run model get.dic iterations to get the DIC.  If \code{parComp = TRUE}, then no DIC is calculated.
#' @param parComp Logical.  If parComp = TRUE, then chains are computed in parallel using the runjags architecture.  This opens n.chains instances of JAGS.  In practice, a max of 4 unique chains can be run due to the way in which JAGS generates initial values.
#' @param autoConv Logical.  If parComp = TRUE and autoConv = TRUE, the chains are run until convergence criteria are met.  Currently, chains are run for a single period equal to \code{burn.in} iterations and monitored for \code{n.iter} iterations.  If PSRF on the three values of pi are lower than 1.05, then the chain is stopped and the chain is run for \code{n.iter} more iterations monitoring all values specified by \code{parameters}.  If the PSRF for any value is higher than 1.05, then the chain is run for \code{burn.in} + \code{n.iter} more iterations and the PSRF is again checked.  This is repeated, at most, \code{max.auto} times.  If PSRF is not met by \code{max.auto} attempts, a warning message is printed and the chains are run \code{n.iter} more times with all parameters monitored.
#' @param max.auto Integer.  Number of subsequent tries to achieve a PSRF of 1.05 on pi.  After \code{max.auto} failures, a warning is thrown and the chain is run \code{n.iter} more times monitoring all specified parameters.
#'
#' @return The function returns a nested list. The first element of the list is a \code{mcmc} object with the samples from the posterior distribution. The second element of the list is a list of summaries (HPD, Mean, etc)
#'
#' @examples
#'
#' model    = 'bl'
#' 
#' ## simulate data
#' ## -------------
#' sim_data = ef_simulateData(n=700, nCov=1, model=model)
#' data     = sim_data$data
#' 
#' ## mcmc parameters
#' ## ---------------
#' mcmc    = list(burn.in=1, n.adapt=10, n.iter=100, n.chains=2)
#' 
#' ## samples
#' ## -------
#' ## help(eforensics)
#' devtools::document(pkg_folder)
#' samples    = eforensics(
#'     w ~ x1.w ,
#'     a ~ x1.a,
#'     mu.iota.m ~ x1.iota.m, 
#'     mu.iota.s ~ x1.iota.s,
#'     ## mu.chi.m  ~ x1.chi.m, 
#'     ## mu.chi.s  ~ x1.chi.s,
#'     data=data,
#'     elegible.voters="N",
#'     model=model, mcmc=mcmc, get.dic=0)
#' 
#' summary(samples)
#' summary(samples, join.chains=T)
#' 
#' @export
## }}}
eforensics   <- function(formula1, formula2, formula3=NULL, formula4=NULL, formula5=NULL, formula6=NULL, data, elegible.voters=NULL, weights=NULL, mcmc, model, parameters="all", na.action="exclude", get.dic = 1000, parComp = T, autoConv = T, max.auto = 10)
{
    ## error handling
    check_mcmc(mcmc)
    options(warn=-1)
    on.exit(options(warn=0))
    ## check if JAGS is installed
    ef_check_jags()

    ## order the formulas (see comment in the order.formulas() function)
    formulas = order.formulas(formula3, formula4, formula5, formula6)
    formula3 = formulas[[1]] %>% as.formula
    formula4 = formulas[[2]] %>% as.formula
    formula5 = formulas[[3]] %>% as.formula
    formula6 = formulas[[4]] %>% as.formula

    ## create a placeholder for weights it if is not provided
    if (is.null(weights)) {data = data %>% dplyr::mutate(weights = 1)}
    
    if(parComp == T){
      eforensics_main_par(formula1, formula2, formula3, formula4, formula5, formula6, data, elegible.voters, weights, mcmc, model, parameters, na.action, get.dic, autoConv, max.auto)
    }else{
      eforensics_main(formula1, formula2, formula3, formula4, formula5, formula6, data, elegible.voters, weights, mcmc, model, parameters, na.action, get.dic)
    }
    
}

eforensics_main   <- function(formula1, formula2, formula3, formula4, formula5, formula6, data, elegible.voters=NULL, weights, mcmc, model, parameters=NULL, na.action="exclude", get.dic = 1000)
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
    if(!is.null(elegible.voters)){
        data = data %>% dplyr::rename(elegible.voters = !!elegible.voters) 
        dat$N = data$elegible.voters
    }else{
        ## check if model use counts, and require elebigle voters
        ## ------------------------------------------------------
        if (stringr::str_detect(model, pattern="bl")) {
            stop("\nThe parameter 'elegible.voters' must be provided for models based on binomial distributions\n\n")
        }
    }
    data = dat
    ## get parameters to monitor
    ## -------------------------
    if(is.null(parameters)) parameters = ef_get_parameters_to_monitor(model)
    if(parameters[1] == 'all') parameters = ef_get_parameters_to_monitor(model, all=TRUE)

    ## get model
    ## ---------
    model.name = model
    model      = get_model(model.name)

    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n','Burn-in: ', mcmc$burn.in, '\n'); cat(msg)
    ## msg <- paste0('\n','Chains: ', mcmc$n.chains, '\n'); cat(msg)
    msg <- paste0('\n','Number of MCMC samples per chain: ', mcmc$n.iter, '\n'); cat(msg)
    msg <- paste0('\n','MCMC in progress ....', '\n'); cat(msg)
    ## ---------------------------------------------------

    ## MCMC
    ## ----
    time.init    = Sys.time()
    cat('\nCompiling the model...\n')     ; sim = rjags::jags.model(file=textConnection(model), data = data, n.adapt=mcmc$n.adapt, n.chain=mcmc$n.chains);
    cat('\nUpdating MCMC (burn-in) ...\n'); stats::update(sim, n.iter = mcmc$burn.in)
    cat('\nDrawing the samples...\n')     ; samples = rjags::coda.samples(model=sim, variable.names=parameters, n.iter=mcmc$n.iter)
    T.mcmc = Sys.time() - time.init
    if(get.dic != 0){
        cat('\nDrawing DIC samples...\n') ; dic.samples = rjags::dic.samples(model=sim, n.iter=get.dic, type = "popt")
    }else{
        dic.samples = NULL
    }
    T.mcmc = Sys.time() - time.init

    if(!is.null(parameters) & "Z" %in% parameters)
        samples = get_Z(samples)
    else
        samples = create_list(samples)
    class(samples) = "eforensics"

    attr(samples, "formula.w") = formula1
    attr(samples, "formula.a") = formula2
    attr(samples, "model")     = model.name
    if (model.name %in% c("rn")) {
        attr(samples, "terms")     = c("alpha", colnames(X.chi.m), colnames(X.iota.m), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
    }else{
        attr(samples, "terms")     = c(colnames(X.chi.m), colnames(X.chi.s), colnames(X.iota.m), colnames(X.iota.s), colnames(Xw), colnames(Xa), "No Fraud", "Incremental Fraud", "Extreme Fraud")
    }
    attr(samples, "dic")       = dic.samples
    ## include alpha in terms name if using rn model

    cat("\n\nEstimation Completed\n\n")
    return(samples)
}

eforensics_main_par   <- function(formula1, formula2, formula3, formula4, formula5, formula6, data, elegible.voters=NULL, weights, mcmc, model, parameters=NULL, na.action="exclude", get.dic = 1000, autoConv = T, max.auto = 10)
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
    if(!is.null(elegible.voters)){
        data = data %>% dplyr::rename(elegible.voters = !!elegible.voters) 
        dat$N = data$elegible.voters
    }else{
        ## check if model use counts, and require elebigle voters
        ## ------------------------------------------------------
        if (stringr::str_detect(model, pattern="bl")) {
            stop("\nThe parameter 'elegible.voters' must be provided for models based on binomial distributions\n\n")
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
    
    ## MCMC
    ## ----
    time.init    = Sys.time()
    ## Run for a prespecified amount of time
    runjags::runjags.options(inits.warning=FALSE, rng.warning=FALSE)
    ## manke sure max.auto is not negative of a fraction. If fraction, get the smaller integer larger than the fraction
    max.auto = max(0, ceiling(max.auto))
    if(autoConv == T){
        trial      = 0
        presim     = runjags::run.jags(model = model, monitor = c("pi"), data = data, n.chains = mcmc$n.chains, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = .5, method = "parallel")
        diagnostic = ef_get_diagnostic(presim)
        while(trial < max.auto & !diagnostic$converged){
            ## Messages about convergence
            ef_print_diagnostic(diagnostic)
            trial      = trial + 1
            cat(paste("\nConvergence diagnostic requirement not met. Extending the chains (attempt ", trial, " of ", max.auto, ") ... \n\n",sep = ""))
            ## extending the chain
            presim = runjags::extend.jags(presim, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = 10, method = "parallel")
            diagnostic = ef_get_diagnostic(presim)
        }
        ## Final sample to keep
        cat("\nBurnin Finished.\nCapturing the samples ...\n\n")
        samples = runjags::extend.jags(presim, add.monitor = c(parameters), burnin = 0, sample = mcmc$n.iter, adapt = mcmc$n.adapt, thin = 1, method = "parallel", jags.refresh = 10)
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
        samples = runjags::run.jags(model = model, monitor = parameters, data = data, n.chains = mcmc$n.chains, burnin = mcmc$burn.in, sample = mcmc$n.iter, adapt = mcmc$n.adapt, jags.refresh = .5, method = "parallel")
        cat("\n\n")
        cat(paste("NOTE: Convergence diagnostics were not computed'. Use results with caution.\n", sep=""))
        cat(paste0("Set 'autoConv=T' to compute diagnostic automatically. See help(eforensics).", sep = ""))
        cat("\n")

    }
    dic.samples = NULL
    T.mcmc = Sys.time() - time.init

    if(!is.null(parameters) & "Z" %in% parameters){
        samples = get_Z(samples[[1]])
    }else{
        samples = create_list(samples[[1]])
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
