#'  Precinct-level data from the 2010 mayoral election in Washington, DC
#'
#' This data set contains precinct-level information from the Washington, DC
#' 2010 mayoral election. Only the key information that is required to
#' estimate fraud in favor of the winner of the elecion (Vincent C. Gray)
#' is provided in the data set. It includes:
#'
#' @format A data frame with 143 rows and 5 variables:
#' \describe{
#'   \item{precinct: }{name of the precinct}
#'   \item{NVoters: }{Number of registed voters}
#'   \item{NValid: }{Number of valid votes}
#'   \item{Votes: }{Total votes for Vincent C. Gray}
#'   \item{a: }{Total abstention}
#' }
#' 
"dc2010"
