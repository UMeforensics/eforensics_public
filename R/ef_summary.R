
#' @export
ef_get_true <- function(sim_data)
{
    true = sim_data$parameters %>%
        tibble::data_frame(Parameter=names(.), True=.)  %>%
        dplyr::mutate(Parameter = paste0(stringr::str_extract(Parameter, 'beta.tau|beta.nu|pi|.*chi.(m|s)|.*iota.(m|s)|n') , '[', stringr::str_extract(Parameter, '[0-9]+') ,']'),
                      Parameter = stringr::str_replace(string=Parameter, pattern="\\[NA\\]", replacement="")) 
    return(true)
}

## {{{ docs }}}

#' Classify observations into fraud distributions
#'
#' This function classifies the data points into one of three distributions: incremental fraud, extreme fraud, and no fraud.
#'
#'
#' @param data a data frame with the data
#' @param samples the output of the function \code{eforensics}  
#'
#'
#' @export

## }}}
ef_classify <- function(data, samples){
    if (!"k.hat" %in% names(samples[[1]])) 
        stop("You need to monitor parameter Z as well in order to classify the data. Use 'parameters=\"all\" in the function eforensics()")

    k.list = list()
    for (i in 1:length(samples))
    {
        k.list[[i]] = samples[[i]]$k.hat
    }
    k.list = k.list %>%
        do.call(cbind,.) %>%
        tibble::as_data_frame(.)  %>% 
        data.table::setnames(., paste0("chain.", 1:length(samples), sep='')) %>%
        base::apply(., 1, function(zi) which.max(c(sum(zi==1), sum(zi==2), sum(zi==3)) ) )
    data = cbind(data, fraud.distribution.idx=k.list) %>%
        dplyr::mutate(fraud.distribution.label = dplyr::case_when(fraud.distribution.idx == 1 ~ "no fraud",
                                                                  fraud.distribution.idx == 2 ~ "incremental fraud",
                                                                  fraud.distribution.idx == 3 ~ "extreme fraud",
                                                                  )) 
    return(data)
}

## {{{ docs }}}
#' Summary 
#'
#' This function return summaries of the posterior distribution estimated by the function \code{eforensics}
#'
#'
#' @param object the output of the function \code{eforensics}
#' @param ... join.chains=TRUE can be used to provide summaries of chains after they are combined together 
#'
#' @export
## }}}
summary.eforensics <- function(object, ...)
{
    args = as.list(match.call())
    if('join.chains' %in% names(args)) {
        join.chains = as.logical(as.character(args$join.chains))
    }else{
        join.chains = FALSE
    }
        

    x = object
    terms = attr(x, "terms")
    if (join.chains) {
        samp = x %>% purrr::map(.x=., ~.x['parameters'][[1]])  %>% do.call(rbind,.)
        HPD = samp %>%
            coda::as.mcmc(.) %>%
            coda::HPDinterval(.) %>%
            data.frame(Parameter = row.names(.), Covariate=terms,., row.names = 1:nrow(.)) %>%
            dplyr::rename(HPD.lower = lower, HPD.upper = upper)  

        samp = samp %>%
            coda::as.mcmc(.) %>%
            summary(.) %>%
            .[[1]] %>%
            data.frame(Parameter = row.names(.), Covariate=terms, ., row.names = 1:nrow(.))  %>%
            dplyr::select(Parameter, Mean, SD) %>%
            dplyr::full_join(., HPD , by=c("Parameter"))  %>%
            dplyr::select(Parameter, Covariate, dplyr::everything()) 
    }else{
        samp = list()
        for (i in 1:length(x))
        {
            HPD = x[[i]]$parameters %>%
                coda::HPDinterval(.) %>%
                data.frame(Parameter = row.names(.), Covariate=terms,., row.names = 1:nrow(.)) %>%
                dplyr::rename(HPD.lower = lower, HPD.upper = upper)  
            tab = x[[i]]$parameters %>%
                summary(.) %>%
                .[[1]] %>%
                data.frame(Parameter = row.names(.), Covariate=terms, ., row.names = 1:nrow(.))  %>%
                dplyr::select(Parameter, Mean, SD) %>%
                dplyr::full_join(., HPD , by=c("Parameter"))  
            samp[[i]] = tab %>%
            dplyr::select(Parameter, Covariate, dplyr::everything()) 
        }
        names(samp) = paste0('Chain ', 1:length(x)) 
    }
    return(samp)
}

#' @export
summary.eforensics_sim_data <- function(object,...)
{
    ## get additional parameters ...
    elegible.voters.provided =FALSE
    args = as.list(match.call())
    if(!'elegible.voters' %in% names(args)) {
        elegible.voters = NULL
    }else{
        elegible.voters = eval(args$elegible.voters)
        elegible.voters.provided =TRUE
    }

    
    sim_data = object
    nrows.data = nrow(sim_data$data)
    par = sim_data$parameters %>%
        as.data.frame(.) %>%
        cbind(., Parameter=row.names(.)) %>%
        dplyr::rename(True=".")  %>%
        tidyr::spread(., key=Parameter, value=True) %>%
        base::replicate(., n=nrows.data, simplify=FALSE) %>%
        base::do.call(base::rbind,.) %>%
        tibble::as_data_frame(.)  %>%
        dplyr::rename_all(dplyr::funs(paste0(., ".True") ))
    sim_data_summ = sim_data$data %>%
        dplyr::bind_cols(., 
                         sim_data$latent %>%
                         base::do.call(base::cbind, .) %>%
                         tibble::as_data_frame(.) %>% 
                         dplyr::rename_all(dplyr::funs(paste0(., ".True") ))
                         ) %>%
        tibble::as_data_frame(.)  %>%
        dplyr::bind_cols(., par) 


    if ('alpha.True' %in% names(sim_data_summ)) {
        sim_data_summ = sim_data_summ %>%
            dplyr::mutate(iota.s.True = iota.m.True^alpha.True,
                          chi.s.True  = chi.m.True^alpha.True) 
    }
    sim_data_summ = sim_data_summ  %>%
        dplyr::mutate(Manufactured.Fraud.True = dplyr::case_when(z.True == 1 ~ 0,
                                                                 z.True == 2 ~ iota.m.True * (1-tau.True),
                                                                     z.True == 3 ~ chi.m.True  * (1-tau.True)) ,
                      Stolen.Fraud.True = dplyr::case_when(z.True == 1 ~ 0,
                                                           z.True == 2 ~ iota.s.True * tau.True * (1-nu.True),
                                                               z.True == 3 ~ chi.s.True * tau.True * (1-nu.True)),
                      Total.Fraud.True = Manufactured.Fraud.True + Stolen.Fraud.True )
    if (elegible.voters.provided) {
        N = sim_data_summ[,N] %>% dplyr::pull(.)
        sim_data_summ = sim_data_summ  %>%
            dplyr::mutate(Manufactured.Fraud.True = N*Manufactured.Fraud.True,
                          Stolen.Fraud.True = N*Stolen.Fraud.True ,
                          Total.Fraud.True = N*Total.Fraud.True)
    }
    return(sim_data_summ)
}

#' @export
print.eforensics <- function(x, ...)
{
    print(summary(x))
    return(summary(x))
}

## {{{ docs }}}
#' Get parameters
#'
#' This function returns the parameters of the model
#'
#'
#' @param samples the output of the function \code{eforensics} 
#'
#' @export
## }}}
ef_get_parameters <- function(samples)
{
    samp = list()
    for (i in 1:length(samples))
    {
        samp[[i]] = samples[[i]]$parameters
    }
    return(samp)   
}


## {{{ docs }}}
#' Compute Fraud Distribution
#'
#' This function returns the estimated posterior distribution of frauds at the individual-level observation
#'
#'
#' @param data the data set used in the estimation
#' @param samples the output of \link{eforensics} function
#' @param model a string with the model used in the estimation (e.g., 'rn', 'bl'). See \link{eforensics} function
#' @param elegible.voters either \code{NULL} (default) or a string with the name of the column containing the number of elegible voters
#' @inheritParams plot_local_fraud 
#' @return It returns a tibble data.frame with the data points classified into no fraud, incremental, or extreme fraud distribution along with the posterior distribution of fraud for each observation
#'
#' @export
## }}}
ef_get_local_fraud <- function(data, samples, model, elegible.voters=NULL, sim_data=NULL)
{
    op.default <- options()
    on.exit(options(op.default), add=TRUE)
    options(warn=-1)

    ## check_local_fraud(samples, model)
    if(!'fraud.distribution.idx' %in% names(data)) dat = ef_classify(data, samples) %>% tibble::as_data_frame(.) 

    ## getting parameters
    ## ------------------
    samples.par            = ef_get_parameters(samples) %>% runjags::combine.mcmc(.) %>% tibble::as_data_frame(.) 
    beta.nu                = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("beta.nu")) 
    beta.tau               = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("beta.tau")) 
    iota.m                 = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.iota.m")) 
    chi.m                  = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.chi.m")) 
    if(model=='bl') {
        iota.s = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.iota.s"))   
        chi.s  = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.chi.s")) 
        sigma.tau = NULL
        sigma.nu  = NULL
    }
    if(model=='rn') {
        alpha     = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("alpha")) %>% dplyr::summarise(mean(alpha)) %>% dplyr::pull(.)
        iota.s    = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.iota.m")) %>% dplyr::mutate(mu.iota.s = mu.iota.m^alpha) %>% dplyr::select(mu.iota.s)   
        chi.s     = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("mu.chi.m"))  %>% dplyr::mutate(mu.chi.s  = mu.chi.m ^alpha) %>% dplyr::select(mu.chi.s)   
        sigma.tau = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("sigma.tau")) 
        sigma.nu  = samples.par  %>% tibble::as_data_frame(.)  %>% dplyr::select(dplyr::contains("sigma.nu")) 
    }

    ## getting data (Xw and Xa) (design matrix)
    ## ------------------------
    getX <- function(formula1, data)
    {
        func.call <- match.call(expand.dots = FALSE)
        mat     = getRegMatrix(func.call, data, weights=NULL, formula_number=1)
        Xw      = mat$X
        return(Xw)
    }
    formula   = attr(samples, 'formula.w')
    Xw = getX(formula, data)
    formula   = attr(samples, 'formula.a')
    Xa = getX(formula, data)

    ## Get fraud distribution
    ## ---------------------------------------------
    Z = dat$fraud.distribution.idx
    N=rep(NA, nrow(dat))
    if (!is.null(elegible.voters) ) {
        if (elegible.voters %in% names(dat) ) {
            N = dat[,elegible.voters]
            N = N %>% dplyr::pull(.)
            cat("\nLocal fraud will be returned in counts.\n")
        }else{
            cat(paste0("\n\nNOTE: The column ", elegible.voters, "is not in the data. The function will return fraud proportions instead of counts.") )
        }
    }else{
        cat("\nLocal fraud will be returned in proportion of votes.\n")
    }
    ## Debug/Monitoring message --------------------------
    msg <- paste0('\n','Computing local fraud distribution.\n\nIt may take some time to compute ...',  '\n'); cat(msg)
    ## ---------------------------------------------------
    doParallel::registerDoParallel(parallel::detectCores() - 1) ## -- Parallel (start) -----------------------------------------
    "%dopar%" <- foreach::"%dopar%"
    dat$fraud.distribution = foreach::foreach(Zi = Z,
                                              Ni = N,
                                              Xwi=split(Xw, 1:nrow(Xw)),
                                              Xai=split(Xa, 1:nrow(Xa)),
                                              .maxcombine = nrow(dat),
                                              .combine = list, .multicombine=TRUE) %dopar% ef_local_fraud_get_distribution(beta.tau, beta.nu, iota.m, iota.s, chi.m, chi.s, Zi, Xwi, Xai, Ni, sigma.tau, sigma.nu, model)
    if (2 %in% unique(dat$fraud.distribution.idx) | 3 %in% unique(dat$fraud.distribution.idx) ) {
        ## Total
        fraud.summary = foreach::foreach(dist = dat$fraud.distribution,
                                         .maxcombine = nrow(dat),
                                         .combine = list, .multicombine=T) %dopar% if(!is.null(dist[1])){c(mean(dist$Total), coda::HPDinterval(coda::as.mcmc(dist$Total)))}else{c(0,0,0)}
        dat = dat %>%
            dplyr::bind_cols(., do.call(rbind,fraud.summary) %>% 
                                tibble::as_data_frame()  %>%
                                dplyr::rename_at(., 1:3, dplyr::funs(c('Total.Fraud.Mean', "Total.Fraud.HPD.lower", "Total.Fraud.HPD.upper")) )
                             ) 
        ## Manufactured
        fraud.summary = foreach::foreach(dist = dat$fraud.distribution,
                                         .maxcombine = nrow(dat),
                                         .combine = list, .multicombine=T) %dopar% if(!is.null(dist[1])){c(mean(dist$Manufactured), coda::HPDinterval(coda::as.mcmc(dist$Manufactured)))}else{c(0,0,0)}
        dat = dat %>%
            dplyr::bind_cols(., do.call(rbind,fraud.summary) %>% 
                                tibble::as_data_frame()  %>%
                                dplyr::rename_at(., 1:3, dplyr::funs(c('Manufactured.Fraud.Mean', "Manufactured.Fraud.HPD.lower", "Manufactured.Fraud.HPD.upper")) )
                             ) 
        ## Stolen
        fraud.summary = foreach::foreach(dist = dat$fraud.distribution,
                                         .maxcombine = nrow(dat),
                                         .combine = list, .multicombine=T) %dopar% if(!is.null(dist[1])){c(mean(dist$Stolen), coda::HPDinterval(coda::as.mcmc(dist$Stolen)))}else{c(0,0,0)}
        dat = dat %>%
            dplyr::bind_cols(., do.call(rbind,fraud.summary) %>% 
                                tibble::as_data_frame()  %>%
                                dplyr::rename_at(., 1:3, dplyr::funs(c('Stolen.Fraud.Mean', "Stolen.Fraud.HPD.lower", "Stolen.Fraud.HPD.upper")) )
                             ) 

    }
    doParallel::stopImplicitCluster()  ## --------------------------- Parallel (stop) -------------------------------------------

    ## get and merge true value if provided
    if (!is.null(sim_data)) {
        dat  = dat %>% dplyr::full_join(., summary(sim_data), by=c()) 
        if (!is.null(elegible.voters) ) {
            if (elegible.voters %in% names(dat)) {
                N = dat[,elegible.voters] %>% dplyr::pull(.)
                dat = dat  %>%
                    dplyr::mutate(Manufactured.Fraud.True = N*Manufactured.Fraud.True,
                                  Stolen.Fraud.True = N*Stolen.Fraud.True ,
                              Total.Fraud.True = N*Total.Fraud.True)
            }
        }
    }

    return(dat)
}
ef_local_fraud_get_distribution <- function(beta.tau, beta.nu, iota.m, iota.s, chi.m, chi.s, Zi, Xwi, Xai, Ni, sigma.tau=NULL, sigma.nu=NULL, model )
{
    if (Zi==1) {
        fraud.density = NULL
    }else{
        mu.tau = as.matrix(beta.tau) %*% matrix(t(Xai))
        mu.nu  = as.matrix(beta.nu ) %*% matrix(t(Xwi))
        if (model=='bl') {
            mu.tau        = 1/(1+exp(-mu.tau))
            mu.nu         = 1/(1+exp(-mu.nu))
        }
        if (model=='rn') {
            ## compute the expectation of a restricted normal distribution
            sigma.tau = sigma.tau %>% dplyr::pull(.)
            a.tau     = (0 - mu.tau)/sigma.tau
            b.tau     = (1 - mu.tau)/sigma.tau
            k.tau     = stats::pnorm(b.tau, 0, 1) - stats::pnorm(a.tau, 0, 1)
            mu.tau    = mu.tau + ( ( stats::dnorm(a.tau) - stats::dnorm(b.tau) ) / k.tau ) * sigma.tau 

            sigma.nu = sigma.nu %>% dplyr::pull(.)
            a.nu     = (0 - mu.nu)/sigma.nu
            b.nu     = (1 - mu.nu)/sigma.nu
            k.nu     = stats::pnorm(b.nu, 0, 1) - stats::pnorm(a.nu, 0, 1)
            mu.nu    = mu.nu + ( ( stats::dnorm(a.nu) - stats::dnorm(b.nu) ) / k.nu ) * sigma.nu 
        }

        if (Zi == 2) {
            manufactured = iota.m * (1-mu.tau)   
            stolen       = iota.s * mu.tau * (1-mu.nu)  
        }
        if(Zi == 3){
            manufactured = chi.m * (1-mu.tau) 
            stolen       = chi.s * mu.tau * (1-mu.nu) 
        }
        if (!is.na(Ni)) {
            manufactured = Ni*manufactured 
            stolen       = Ni*stolen 
        }
        total         = manufactured + stolen
        fraud.density = data.frame(manufactured=manufactured, stolen=stolen, total=total)  %>% tibble::as_data_frame(.) 
        names(fraud.density) = c('Manufactured', 'Stolen', 'Total')
    }
    return(fraud.density)
}

## =====================================================
## plots
## =====================================================
parameter.rename  <- function(parameter)
{
    parameter = parameter  %>%  
    ## stringr::str_replace(string=., pattern="\\[", replacement="_[") %>%
    ## stringr::str_replace(string=., pattern="]", replacement="}") %>%
    stringr::str_replace(string=., pattern=".tau\\[", replacement="[tau*") %>% 
    stringr::str_replace(string=., pattern=".nu\\[", replacement="[nu*")  %>% 
    stringr::str_replace(string=., pattern=".chi.", replacement="[chi]^")  %>% 
    stringr::str_replace(string=., pattern=".iota.", replacement="[iota]^") %>%
    stringr::str_replace_all(string=., pattern=" ", replacement="~") 
    ## stringr::str_replace(string=., pattern="\\.w~", replacement="[w]~") %>% 
    ## stringr::str_replace(string=., pattern="\\.a~", replacement="[a]~") 
    return(parameter)
}
## {{{ docs }}}
#' Plot summary of posterior distribution
#'
#' Plot summary of the samples from posterior distribution produced by \code{\link{eforensics}} function
#'
#'
#' @inheritParams summary.eforensics
#' @param true used when data was generated by simulation produced by \code{\link{ef_simulateData}} function. The value must the be output of the function \code{\link{ef_get_true()}}
#' @param parse boolean, if \code{TRUE}, text in the tick marks of the y-axis is parsed.
#' @param plots a string vector with any combination of "Abstention and Vote", "Local Fraud Probabilities", and "Fraud Probability". It defines with plot will be included in the figure. Default includes all.
#' @param title a string with the title of the plot
#' @param subtitle a string wiht the subtitle of the plot
#' @param xlab a string with the lable of the x-axis
#' @param ylab a string with the lable of the y-axis
#'
#'
#' @export
## }}}
ef_plot <- function(samples, true=NULL, parse=FALSE, plots=c("Abstention and Vote", "Local Fraud Probabilities", "Fraud Probability" ), title=NULL, subtitle=NULL, xlab=NULL, ylab=NULL)
{
    options(warn=-1)
    on.exit(options(warn=0))
    
    ## creating summary table
    tab = summary(samples, join.chains=T)
    if (!is.null(true)) {
        ## joining true if 'True' value of the parameters are provided
        tab = tab %>% dplyr::left_join(., True)    
    }

    tab = tab %>%
        dplyr::filter( ! (stringr::str_detect(Parameter, pattern="beta.*[1]") & stringr::str_detect(Covariate, pattern="Intercept"))  ) %>% 
        dplyr::mutate(Group = dplyr::case_when(stringr::str_detect(Parameter, pattern="pi")~"Fraud Probability",
                                               stringr::str_detect(Parameter, pattern="beta.(nu|tau)")~"Coefficients: Abstention and Vote for the winner",
                                               TRUE~"Coefficients: Local Fraud Probabilities"
                                               ),
                      Parameter = dplyr::case_when(stringr::str_detect(Parameter, pattern="nu") ~ "votes to the winner",
                                                   stringr::str_detect(Parameter, pattern="tau") ~ "abstension",
                                                   stringr::str_detect(Parameter, pattern="iota.s") ~ "Incremental/Stolen",
                                                   stringr::str_detect(Parameter, pattern="iota.m") ~ "Incremental/Manufactured",
                                                   stringr::str_detect(Parameter, pattern="chi.s") ~ "Extreme/Stolen",
                                                   stringr::str_detect(Parameter, pattern="chi.m") ~ "Extreme/Manufactured",
                                                   ## stringr::str_detect(Parameter, pattern="iota.s") ~ "Incremental Fraud/Stolen votes",
                                                   ## stringr::str_detect(Parameter, pattern="iota.m") ~ "Incremental Fraud/Manufactured votes",
                                                   ## stringr::str_detect(Parameter, pattern="chi.s") ~ "Extreme Fraud/Stolen votes",
                                                   ## stringr::str_detect(Parameter, pattern="chi.m") ~ "Extreme Fraud/Manufactured votes",
                                                   ## stringr::str_detect(Parameter, pattern="pi.1.") ~ "expression(p[1])",
                                                   ## stringr::str_detect(Parameter, pattern="pi.2.") ~ "expression(p[2])",
                                                   ## stringr::str_detect(Parameter, pattern="pi.3.") ~ "expression(p[3])",
                                                   TRUE ~ Parameter %>% as.character),
                      label = paste0(Covariate, " (", Parameter, ")"),
                      label = stringr::str_replace(string=label, pattern=".Intercept.", replacement="Intercept"),
                      )  %>% 
        dplyr::filter(stringr::str_detect(Group, pattern=paste0(plots, collapse="|") )) 
    if (parse) {
        tab = tab %>% dplyr::mutate(label =  parameter.rename(label))
    }else{
        tab = tab %>% dplyr::mutate(label = stringr::str_replace_all(string=label, pattern="\\(pi.[0-9].\\)", replacement="") ) 
    }
    ## plot
    ## ----
    if (is.null(xlab )) {xlab = "Posterior Expectation"}
    if (is.null(ylab )) {ylab = ""}
    g = tab %>%
        ggplot2::ggplot(.)  +
        ggplot2::geom_point(ggplot2::aes(x=Mean, y=label, shape="Posterior Mean", colour='Posterior Mean'), size=2, alpha=1)   +
        ggplot2::geom_errorbarh(ggplot2::aes(y = label, xmin=HPD.lower, xmax=HPD.upper, linetype="95% HPD interval"), height=.1) +
        ggplot2::xlab(xlab)   +
        ggplot2::ylab(ylab)   +
        ggplot2::scale_linetype(name="")  +
        ggplot2::scale_shape(name="")   +
        ggplot2::facet_wrap( ~ Group, ncol = 1, scales='free' )  
    if (!is.null(true)) {
        g = g +
            ggplot2::geom_point(ggplot2::aes(x=True, y=label, shape="True", colour="True"), size=2, alpha=.6)    +
            ggplot2::scale_colour_manual(values = c("Posterior Mean" = "black", "True"='red'), name="")  
    }else{
        
        g = g +
            ggplot2::scale_colour_manual(values = c("Posterior Mean" = "black"), name="")  
    }
    if (parse) {
        g = g + 
            ggplot2::scale_y_discrete(labels=function(l)  parse(text=l)) 
    }
    if (!is.null(title)) {
        g = g +
            ggplot2::ggtitle(label=title, subtitle="") 
    }
    if (!is.null(title) & !is.null(subtitle)) {
        g = g +
            ggplot2::ggtitle(label=title, subtitle=subtitle) 
    }
    ## theme
    g = g +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top") + 
        ggplot2::theme(strip.background = ggplot2::element_rect(colour="white", fill="white"),
                       strip.text.x = ggplot2::element_text(size=11, face='bold', hjust=0),
                       strip.text.y = ggplot2::element_text(size=11, face="bold", vjust=0)) 
    return(g)
}

## put in the package
## ------------------
## {{{ docs }}}

#' Plot samples from eforensics
#'
#' This function creates a plot with posterior distribution of fraud in a given election unit (ballot box, pooling place, etc)
#'
#' @param data either the data set used to estimate fraud or the output of the function \link{ef_get_local_fraud}
#' @param samples eigher \code{NULL} (default) (which requires \code{data} to be the output of the function \link{ef_get_local_fraud}) or the outuput of the function \link{eforensics}
#' @param model either \code{NULL} (default) (which requires \code{data} to be the output of the function \link{ef_get_local_fraud}) or a string with the name of the model used to estimate fraud (see \link{eforensics})
#' @param row an integer with the row number of data to plot
#' @param by.types boolean, if \code{TRUE} incremental, extreme, and total fraud are displayed. Otherwise, only distribution of total fraud is shown
#' @param election.unit a string with the name of the column that contains a label that identifies the election unit. Default \code{NULL}
#' @param plot.mean boolean, if \code{TRUE} the posterior average of the distribution is also displayed as well as the 95\% HPD
#' @param legend.position a string with the position of the legend when posterior mean is displayed. Accepts \code{top}, \code{bottom}, \code{left}, \code{right}
#' @param sim_data the output of the function \link{ef_simulateData}. It is used only if the fraud was estimated using simulated data. If provided, the plot also display the true fraud value. Default \code{NULL}
#' @param title a string with the title of the plot. Default \code{NULL}
#' @param subtitle a string with the subtitle of the plot. Default \code{NULL}
#'
#' @export

## }}}
plot_local_fraud <- function(data, samples=NULL, model=NULL, row=NULL, election.unit=NULL, title=NULL, subtitle=TRUE, plot.mean=TRUE, legend.position='bottom', sim_data=NULL, by.types=TRUE)
{
    ## check if data provided contains classification and if not, classify it
    ## ----------------------------------------------------------------------
    if (!("fraud.distribution.idx" %in% names(data) & "fraud.distribution.label" %in% names(data))) {
        if (is.null(samples) | is.null(model)) {
            stop("\n\n The 'data' provided does not contain classification of observations. Columns 'fraud.distribution.idx' and 'fraud.distribution.label' must be in the data. Either provide 'samples' and 'model' or use the function ef_get_local_fraud() to classify the data first.\n\n")
        }
        data = ef_get_local_fraud(data, samples, model)
    }
    ## checke if row was provided and if it was, if that case was faudulent
    ## --------------------------------------------------------------------
    if (is.null(row)) {
        stop("\n\nThe row number with the observation to plot must be provided. Use the parameter 'row'\n\n")
    }else{
        if (data$fraud.distribution.idx[row]==1) return(cat("\n\nObservation in row ", row, " in the data was classified as a non-fraudulent case.\n\n") )
    }

    ## if simulated data is provided, get summary to plot true values
    ## --------------------------------------------------------------
    if (!is.null(sim_data)) {
        ## Debug/Monitoring message --------------------------
        msg <- paste0('\n','Getting true parameter values ...',  '\n'); cat(msg)
        ## ---------------------------------------------------
        true = summary(sim_data)
        data = data %>%
            dplyr::left_join(., true) 
    }
    

    data = data %>% dplyr::filter(dplyr::row_number()==row)
    if (by.types) {
        g = plot_frauds(data, plot.mean, sim_data)
    }else{
        g = plot_frauds_total(data, plot.mean, sim_data)
    }
    g = g +  
        ggplot2::theme(legend.position = legend.position) 



    ## plot title and subtitle
    if (!is.null(title)) {
        g = g + ggplot2::ggtitle(title)
    }else{
        title=''
    }
    if (subtitle) {
        subtitle = paste0("(", data$fraud.distribution.label, ")") 
        g = g + ggplot2::ggtitle(title,subtitle=subtitle)
    }
    if (!is.null(election.unit)) {
        if (election.unit %in% names(data)) {
            election.unit = paste0(" Election unit: ", data[,election.unit] %>% dplyr::pull(.)) 
            g = g + ggplot2::annotate("text", label = election.unit, x = -Inf, y = Inf, hjust=0, vjust=1.5 )
        }
    }
    return(g)
}
plot_frauds_total <- function(tab, plot.mean, sim_data)
{
    g =  tab  %>%
        dplyr::select(fraud.distribution)  %>%
        dplyr::pull(.) %>%
        .[[1]] %>% 
        dplyr::select(Total)  %>% 
        ggplot2::ggplot(.) +
        ggplot2::geom_histogram(ggplot2::aes(x=Total, y=..density..),fill="lightsteelblue1", colour='white') +
        ggplot2::geom_density(ggplot2::aes(x=Total), fill="#00000044", adjust=1, alpha=.2, colour="white") +
        ggplot2::theme_bw()+
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) 
    if (plot.mean & !is.null(sim_data)) {
        g = g +
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.Mean,  linetype="Mean", colour='Mean'))+
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.HPD.lower,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.HPD.upper,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.True,  linetype="True", colour='True'))+
            ggplot2::scale_linetype_manual(values=c( "95% HPD"='dashed', "Mean"='solid', 'True'='solid'), name='') +
            ggplot2::scale_colour_manual(values = c( "95% HPD"='red', "Mean"='red', 'True'='black'), name='') 
    }
    if (plot.mean & is.null(sim_data)) {
        g = g +
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.Mean,  linetype="Mean", colour='Mean'))+
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.HPD.lower,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab, ggplot2::aes(xintercept=Total.Fraud.HPD.upper,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::scale_linetype_manual(values=c( "95% HPD"='dashed', "Mean"='solid'), name='') +
            ggplot2::scale_colour_manual(values = c( "95% HPD"='red', "Mean"='red'), name='') 
    }
    return(g)
}
plot_frauds <- function(tab, plot.mean, sim_data)
{
    g = tab %>% 
        ## dplyr::filter(dplyr::row_number()==row.number)   %>%
        dplyr::select(fraud.distribution)  %>%
        dplyr::pull(.) %>%
        .[[1]] %>% 
        tidyr::gather(key = Fraud, value=value)  %>% 
        ggplot2::ggplot(.) +
        ggplot2::geom_histogram(ggplot2::aes(x=value, y=..density..),fill="lightsteelblue1", colour='white') +
        ggplot2::geom_density(ggplot2::aes(x=value), fill="#00000044", adjust=1, alpha=.2, colour='white') +
        ggplot2::theme_bw()+
        ggplot2::facet_wrap( ~ Fraud, ncol = 1, scales='free',labeller=ggplot2::label_parsed)  +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_fill_brewer(palette='Blues', name="Fraud")+
        ggplot2::theme(legend.position = "top") +
        ggplot2::theme(strip.background = ggplot2::element_rect(colour="white", fill="white"),
                       strip.text.x = ggplot2::element_text(size=12, face='bold', hjust=0),
                       strip.text.y = ggplot2::element_text(size=12, face="bold", vjust=0)) 
    if (plot.mean & !is.null(sim_data)) {
        tab2 = tab %>%
            dplyr::select(dplyr::contains("Manufactured"), dplyr::contains("Stolen"), dplyr::contains("Total")) %>%
            tidyr::gather(key = Fraud, value=value) %>%
            tidyr::separate(., col=Fraud, into=c("Fraud", "stat"), sep=".Fraud.") %>%
            tidyr::spread(., key=stat, value=value)
        g = g +
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=Mean,  linetype="Mean", colour='Mean'))+
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=HPD.lower,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=HPD.upper,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=True,  linetype="True", colour='True'))+
            ggplot2::scale_linetype_manual(values=c( "95% HPD"='dashed', "Mean"='solid', 'True'='solid'), name='') +
            ggplot2::scale_colour_manual(values = c( "95% HPD"='red', "Mean"='red', 'True'='black'), name='')
    }
    if (plot.mean & is.null(sim_data)) {
        tab2 = tab %>%
            dplyr::select(dplyr::contains("Manufactured"), dplyr::contains("Stolen"), dplyr::contains("Total")) %>%
            tidyr::gather(key = Fraud, value=value) %>%
            tidyr::separate(., col=Fraud, into=c("Fraud", "stat"), sep=".Fraud.") %>%
            tidyr::spread(., key=stat, value=value)
        g = g +
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=Mean,  linetype="Mean", colour='Mean'))+
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=HPD.lower,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::geom_vline(data=tab2, ggplot2::aes(xintercept=HPD.upper,  linetype="95% HPD", colour="95% HPD"))+
            ggplot2::scale_linetype_manual(values=c( "95% HPD"='dashed', "Mean"='solid'), name='') +
            ggplot2::scale_colour_manual(values = c( "95% HPD"='red', "Mean"='red'), name='')
    }
    return(g)
}

