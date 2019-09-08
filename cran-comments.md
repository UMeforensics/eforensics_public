## Test environments
* local ubuntu 18.04  R 3.4.4
* OS X install (on travis-ci) R 3.4.4
* win-builder (on travis-ci, devel and release) R 3.4.4

## devtools::check() ubuntu results

0 errors | 0 warnings | 0 note

## R CMD check ubuntu results

0 errors | 0 warnings | 1 note

NOTE refers to the spurious complaint regarding other packages:

* checking dependencies in R code ... NOTE
Namespaces in Imports field not imported from:
  ‘LCA’ ‘rjags’
  All declared Imports should be used.

* This is a new release.
