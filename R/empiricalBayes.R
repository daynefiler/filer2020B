#' @name empiricalBayesFuncs
#' @title Functions for performing empirical Bayesian calculations
#' @description Collection of functions for performing empirical Bayesian
#' calculations
#' @inheritParams pkgParams
#' @param empPMAR Named numeric, empirical estimates of expected PMAR values
#' @param empGenoFreq Named numeric, empirical estimates of genotype pair
#' frequencies
#' @param empGeno Character, empirical genotype calls
#' @param ffPrior Numeric of length 1, fetal fraction prior
#' @param ffPriorWt Numeric of length 1, the weight given to the fetal fraction
#' prior
#' @param qPrior Numeric of length 1, average minor allele frequency prior
#' @param qPriorWt Numeric of length 1, the weight given to the q prior
#' @param pltTrain Logical of length 1, plot genotypes at each training
#' iteration when TRUE
#' @details
#' ffPriorWt can be thought of as the number of observations supportint the
#' fetal fraction prior
NULL

#' @describeIn empiricalBayesFuncs Update the emperic genotype call
#' @importFrom stats dbinom
#' @importFrom dlfUtils l1int logSumExp
#' @export

cfUpdateEmpGeno <- function(mar, dep, empPMAR, empGenoFreq) {
  l1int(type.convert(mar))
  l1int(type.convert(dep))
  data("GenoMeta", package = "filer2020B", envir = environment())
  stopifnot(all(hasName(empPMAR, GenoMeta$name)))
  stopifnot(all(hasName(empGenoFreq, GenoMeta$name)))
  empPMAR <- empPMAR[GenoMeta$name]
  empGenoFreq <- empGenoFreq[GenoMeta$name]
  ll <- dbinom(mar, dep, empPMAR, log = TRUE) + log(empGenoFreq)
  pp <- logSumExp(ll)
  pp[which.max(pp)]
}

#' @describeIn empiricalBayesFuncs Update the emperic PMAR values
#' @export

cfUpdateEmpPMAR <- function(marVec, depVec, empGeno, ffPrior, ffPriorWt) {
  data("GenoMeta", package = "filer2020B", envir = environment())
  a <- ffPriorWt*cfExpPMAR(ffPrior)[GenoMeta$name]
  b <- ffPriorWt - a
  marSum <- tapply(marVec, empGeno, sum)[GenoMeta$name]
  depSum <- tapply(depVec, empGeno, sum)[GenoMeta$name]
  marSum[is.na(marSum)] <- 0; names(marSum) <- GenoMeta$name
  depSum[is.na(depSum)] <- 0; names(depSum) <- GenoMeta$name
  (marSum + a)/(depSum + a + b)
}

#' @describeIn empiricalBayesFuncs Update the emperic genotype frequencies
#' @export

cfUpdateEmpGenoFreq <- function(empGeno, qPrior, qPriorWt) {
  data("GenoMeta", package = "filer2020B", envir = environment())
  d <- cfExpGenoFreq(qPrior)[GenoMeta$name]*qPriorWt
  gtN <- tapply(rep(1, length(empGeno)), empGeno, sum)[GenoMeta$name]
  gtN[is.na(gtN)] <- 0; names(gtN) <- GenoMeta$name
  (d + gtN)/sum(d + gtN)
}

#' @describeIn empiricalBayesFuncs Make empiric genotype calls
#' @importFrom dlfUtils l1log l1num l1int
#' @export

cfEmpGeno <- function(marVec, depVec, qPrior, ffPrior, nchng = 5, its = 50,
                      ffPriorWt = length(marVec)/5, qPriorWt = 1000,
                      pltTrain = FALSE) {

  stopifnot(is.integer(type.convert(marVec)))
  stopifnot(is.integer(type.convert(depVec)))
  l1int(type.convert(nchng))
  l1int(type.convert(its))
  l1num(ffPrior)
  l1num(ffPriorWt)
  l1num(qPrior)
  l1num(qPriorWt)
  l1log(pltTrain)

  ev <- list(empGenoFreq = cfExpGenoFreq(qPrior), empPMAR = cfExpPMAR(ffPrior))
  i <- 0L
  empGeno1 <- rep("null", length(marVec))
  repeat {
    empGeno0 <- empGeno1
    empGeno <- mapply(cfUpdateEmpGeno, mar = marVec, dep = depVec,
                      MoreArgs = ev, SIMPLIFY = TRUE)
    empGeno1 <- names(empGeno)
    if (pltTrain) {
      cfPltFreqHist(marVec, depVec, empGeno1, legend = TRUE)
      title(ylab = sprintf("Iteration: %d", i), line = 0)
    }
    ev <- list(empGenoFreq   = cfUpdateEmpGenoFreq(empGeno1, qPrior, qPriorWt),
               empPMAR = cfUpdateEmpPMAR(marVec, depVec, empGeno,
                                         ffPrior, ffPriorWt))
    if (sum(empGeno0 != empGeno1) < nchng || i == its) break
    i <- i + 1L
  }

  empGeno

}




