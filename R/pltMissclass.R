#' @import graphics
#' @import grDevices
#' @export

pltMissclass <- function() {
  weitzman <- function(ff, dep) {
    d <- lapply(c(0.5, 0.5 + ff/2), function(x) dbinom(0:dep, dep, prob = x))
    sum(pmin(d[[1]], d[[2]]))
  }

  deps <- seq(50, 10000, 50)
  ff <- seq(0.05, 0.25, 0.0125)
  ffN <- length(ff)
  res <- lapply(ff, function(x) sapply(deps, weitzman, ff = x))
  data(GenoMeta, envir = environment())
  blues <- colorRampPalette(GenoMeta$color[c(5, 1)])(length(res))

  par(mar = c(4, 4, 1, 1) + 0.1)
  plot.new()
  ofig <- par('fig')
  plot.window(range(deps), ylim = 0:1)
  # abline(h = seq(0.05, 0.95, 0.05), lty = "dashed", col = "lightgray")
  abline(h = 0.05, lty = "dashed", col = "darkgray", lwd = 1)
  mapply(points,
         y = res,
         col = blues,
         MoreArgs = list(x = deps, type = "l", lwd = 1))
  axis(side = 1)
  axis(side = 2)
  title(xlab = "Sequencing depth", ylab = "Error rate (ABab vs. ABbb)")

  par(mar = c(3, 0, 0, 0),
      new = TRUE,
      fig = c(0.5, 0.9, 0.7, 0.9),
      mgp = c(2, 1, 0),
      cex = 0.75)
  plot.new()
  plot.window(xlim = c(0, ffN), ylim = c(0, 1))
  rect(0:(ffN - 1), 0, 1:(ffN), 1, border = NA, col = blues)
  axis(side = 1,
       at = seq(0.5, 16.5, 1),
       labels = ifelse(ff %% 0.05 == 0, sprintf("%0.2f", ff), NA),
       lwd = 0,
       lwd.ticks = 1)
  title(xlab = "Fetal fraction")

  par(fig = ofig, new = TRUE, cex = 1)

}
