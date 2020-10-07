#' @name cfPltPMARHist
#' @title Plot histograms of PMAR values by genotype
#' @description Plot histograms of PMAR values by genotype
#' @details
#' When ff is not NULL, add normal approximation lines.
#' @inheritParams pkgParams
#' @param legend logical of length 1, add a genotype legend when TRUE
#' @param ff Numeric of length 1, fetal fraction
#' @importFrom dlfUtils allSameLength col2alpha l1log line2user l1num
#' @importFrom graphics title hist par legend
#' @export

cfPltFreqHist <- function(marVec, depVec, gtpVec, legend = FALSE, ff = NULL) {
  stopifnot(allSameLength(marVec, depVec, gtpVec))
  stopifnot(is.integer(type.convert(marVec)))
  stopifnot(is.integer(type.convert(depVec)))
  stopifnot(l1log(legend))
  stopifnot(is.null(ff) || l1num(ff))
  data("GenoMeta", package = "filer2020B", envir = environment())
  stopifnot(is.character(gtpVec) || all(gtpVec %in% GenoMeta$name))
  f <- marVec/depVec
  b <- seq(0, 1, 0.005)

  if (!is.null(ff)) {
    expPmar <- cfExpPMAR(ff)[-c(1, 7)]
    expPmarSd <- cfExpSdPMAR(ff, dep = median(depVec))[-c(1, 7)]
    mult <- table(gtpVec)[GenoMeta$name[-c(1, 7)]]
    expCount <- function(mn, sd, mult) {
      diff(pnorm(b, mean = mn, sd = sd))*mult
    }
    counts <- mapply(expCount, expPmar, expPmarSd, mult)
    ylim <- range(counts)
  } else {
    ylim <- NULL
  }

  par(oma = c(ifelse(legend, 3, 0), 0, 0, 0), mar = c(4, 2, 1, 1))
  addHist <- function(x, gVec, g, col) {
    hist(x[gVec == g], breaks = b, col = col, add = TRUE, border = NA)
  }
  cols <- sapply(GenoMeta$color, col2alpha)
  xcld <- GenoMeta$name[c(1, 7)]
  hist(f[!gtpVec %in% xcld], breaks = b, border = NA,
       ann = FALSE, yaxt = "n", ylim = ylim)
  title(xlab = "Observed PMAR")
  mapply(addHist, g = GenoMeta$name, col = cols,
         MoreArgs = list(x = f, gVec = gtpVec))

  if (!is.null(ff)) {
    for (i in seq(5)) {
      tx <- (b[-1] + 0.005)[counts[ , i] > 0.5]
      ty <- counts[counts[ , i] > 0.5, i]
      lines(x = tx, y = ty, col = GenoMeta$color[i + 1])
    }
  }

  if (legend) {
    legend(x = 0.5,
           y = line2user(line = 1, side = 1, outer = TRUE),
           legend = GenoMeta$name,
           fill = GenoMeta$color,
           horiz = TRUE,
           xpd = NA,
           bty = "n",
           xjust = 0.5)
  }

}





