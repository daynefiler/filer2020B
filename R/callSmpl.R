#' @export

callSmpl <- function(a, d, N, f) {
  empGeno <- cfEmpGeno(marVec = a,
                       depVec = d,
                       qPrior = 0.005,
                       qPriorWt = N*2,
                       ffPrior = f,
                       ffPriorWt = N*2)
  empFF <- cfEstFF(marVec = a, depVec = d, gtpVec = names(empGeno))
  calls <- mapply(cfCallGeno, mar = a, dep = d, MoreArgs = list(ff = empFF))
  list(empFF, names(calls), calls)
}
