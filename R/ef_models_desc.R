## {{{ docs }}}

#' Election Forensics Finite Mixture Model Specifications
#'
#' There are currently 3 supported models included in this package.  Each of these strings can be used as the \code{model} argument in \code{\link{eforensics}}
#'
#'
#' @param qbl Estimates the quasi-binomial logistic model of frauds.  Requires the count of eligible voters for each observation (\code{eligible voters != NULL}).  Models the number of abstentions and number of votes for the winner using a quasi-binomial likelihood.  This specification helps to account for overdispersion (as compared to a binomial likelihood) in the counts of abstentions and votes for the winner in both fraudulent and non-fraudulent settings.
#'  
#' Under \code{qbl}, the mixing parameters and coefficients on covariates are returned.  This model also returns estimates for the proportion of manufactured and stolen votes for each observation.  The posterior mean, 95 percent HPD, and quantiles of manufactured votes for each observation can be attained using \code{attr(foo, "frauds")$Manufactured}.  The posterior mean, 95 percent HPD, and quantiles of stolen votes for each observation can be attained using \code{attr(foo, "frauds")$Stolen}.  This model is the default and is the recommended specification for best estimating the model of election frauds.
#' 
#' @param bl Estimates the binomial logistic model of frauds.  Requires the count of eligible voters for each observation (\code{eligible voters != NULL}).  Models the number of abstentions and number of votes for the winner using a binomial likelihood.  Unlike \code{qbl}, the number of abstentions and votes for the winner is assumed to follow a strict binomial specification with corresponding mean and variance.  Under overdispersion in counts, this model can lead to overestimation of the proportion of observations that are fraudulent.
#'   
#' Under \code{bl}, the mixing parameters and coefficients on covariates are returned.  This model also returns estimates for the proportion of manufactured and stolen votes for each observation.  The posterior mean, 95 percent HPD, and quantiles of manufactured votes for each observation can be attained using \code{attr(foo, "frauds")$Manufactured}.  The posterior mean, 95 percent HPD, and quantiles of stolen votes for each observation can be attained using \code{attr(foo, "frauds")$Stolen}.
#' 
#' @param rn Estimates the restricted normal model of frauds.
#' 
#' @usage NULL   
#'
#' @export
ef_models_desc <- function(){}

