# required (by devtools) to link the cpp code 
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%

.onAttach <- function(libname, pkgname){
    packageStartupMessage('

 ---------------------------------------
 Election Forensics Package (eforensics)
 ---------------------------------------

 Authors:

 Supported by NSF grant SES 1523355


 ')
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "contrasts",
                                                        "lower", "upper",
                                                        "Mean", "SD", "Parameter",
                                                        "dist",
                                                        "model ",
                                                        "mu.chi.m",
                                                        "mu.chi.s",
                                                        "mu.iota.m",
                                                        "mu.iota.s",
                                                        "Ni",
                                                        "samples.par",
                                                        "Xai",
                                                        "Xwi",
                                                        "Zi",
                                                        "model",
                                                        "fraud.distribution",
                                                        "value",
                                                        "Fraud",
                                                        "..density..",
                                                        "label_parsed",
                                                        "chi.m.True",
                                                        "iota.m.True",
                                                        "alpha.True",
                                                        "Manufactured.Fraud.True",
                                                        "Stolen.Fraud.True",
                                                        "HPD.lower",
                                                        "HPD.upper",
                                                        "stat",
                                                        "Total",
                                                        "Total.Fraud.True",
                                                        "Total.Fraud.Mean",
                                                        "Total.Fraud.HPD.lower",
                                                        "Total.Fraud.HPD.upper",
                                                        "True"
                                                        ))
