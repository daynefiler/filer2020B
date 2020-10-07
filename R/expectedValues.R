#' @name expectedValues
#' @title Calculate expected values
#' @description Collection of functions for calculating various expected values
#' @inheritParams pkgParams
#' @details
#' Default seqError of 0.003 taken from PMID30026539 (mean + 1SD to over-
#' estimate).
NULL

#' @describeIn expectedValues Expected proportion of minor allele reads
#' @export

cfExpPMAR <- function(ffVec, seqError = 0.003, clpsList = TRUE) {
  n <- length(ffVec)
  l <- list(rep(0 + seqError, n),
            ffVec/2,
            (1 - ffVec)/2,
            rep(1/2, n),
            (1 + ffVec)/2,
            1 - ffVec/2,
            rep(1 - seqError, n))
  data("GenoMeta", package = "filer2020B", envir = environment())
  names(l) <- GenoMeta$name
  if (n == 1 && clpsList) l <- unlist(l)
  l
}

#' @describeIn expectedValues Expected standard deviation on proportion of minor
#' allele reads
#' @export

cfExpSdPMAR <- function(ffVec, dep, seqError = 0.003, clpsList = TRUE) {
  n <- length(ffVec)
  binSD <- function(x, n) sqrt(x*(1 - x)/n)
  l <- lapply(cfExpPMAR(ffVec, seqError = seqError), binSD, n = dep)
  if (n == 1 && clpsList) l <- unlist(l)
  l
}

#' @describeIn expectedValues Expected proportion of sites with unique fetal
#' alleles
#' @export

cfExpUnqHet <- function(qVec) {
  pVec <- 1 - qVec
  pVec^2 - pVec^3 + qVec^2 - qVec^3
}

#' @describeIn expectedValues Expected proportion of genotype pairs
#' @export

cfExpGenoFreq <- function(qVec, clpsList = TRUE) {
  n <- length(qVec)
  pVec <- 1 - qVec
  l <- list(pVec^3,
            pVec^2 - pVec^3,
            pVec^2 - pVec^3,
            pVec - pVec^2,
            qVec^2 - qVec^3,
            qVec^2 - qVec^3,
            qVec^3)
  data("GenoMeta", package = "filer2020B", envir = environment())
  names(l) <- GenoMeta$name
  if (n == 1 && clpsList) l <- unlist(l)
  l
}





