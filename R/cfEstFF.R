#' @title Estimate fetal fraction
#' @description Estimate fetal fraction
#' @inheritParams pkgParams
#' @export

cfEstFF <- function(marVec, depVec, gtpVec) {
  ## Need to do the test -- is it better to average proportions or calculate
  ## one proportion with all the data... similar to an average weighted on the
  ## depth at each unique allele
  pmar <- marVec/depVec
  ff1 <- pmar[gtpVec == "AAab"]
  ff2 <- 1 - pmar[gtpVec == "BBab"]
  ff <- median(c(ff1, ff2))*2
  ff
}