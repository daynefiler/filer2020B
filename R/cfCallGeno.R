#' @name cfCallGeno
#' @title Calculate maximum likelihood genotype
#' @description Calculate maximum likelihood genotype
#' @details 
#' Add later.
#' @inheritParams pkgParams
#' @return Posterior probability for the maximum likelihood genotype
#' @importFrom stats dbinom
#' @importFrom dlfUtils logSumExp l1num
#' @export

cfCallGeno <- function(mar, dep, ff) {
  stopifnot(l1int(type.convert(mar)))
  stopifnot(l1int(type.convert(dep)))
  stopifnot(l1num(ff))
  genos <- cfExpPMAR(ff)
  ll <- dbinom(mar, dep, genos, log = TRUE)
  pp <- logSumExp(ll)
  pp[which.max(pp)]
}
