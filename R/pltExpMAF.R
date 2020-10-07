#' @import graphics
#' @export

pltExpMAF <- function(dep) {

  data(GenoMeta, envir = environment())
  rng <- c(0.025, 0.25)
  f <- seq(rng[1], rng[2], 0.005)
  expPmar <- cfExpPMAR(f)

  par(mar = c(4, 4, 1, 1) + 0.1)
  plot.new()
  plot.window(xlim = c(rng[1], rng[2]), ylim = c(0, 1))
  abline(v = seq(0, 0.5, 0.05), lty = "dashed", col = "gray80")
  for(i in seq_along(expPmar)) {
    PmarSd <- cfExpSdPMAR(f, dep)
    polygon(x = c(f, rev(f)),
            y = c(expPmar[[i]] + 1.96*PmarSd[[i]],
                  rev(expPmar[[i]] - 1.96*PmarSd[[i]])),
            col = col2alpha(GenoMeta$color[i], 0.25),
            border = NA)
  }
  axis(side = 2); # mtext(side = 2, line = 3, "95% CI on Expected MAF")
  axis(side = 1,
       at = c(9.02, 14.75, 19.58)/100,
       labels = c("10-15", "25-30", "30-40"),
       tcl = -0.8,
       mgp = c(3, 2, 0),
       col = "gray50",
       col.axis = "gray50",
       xpd = NA)
  axis(side = 1, tcl = -0.4)
  title(ylab = "PMAR 95% CI",
        xlab = "Fetal fraction & gestational age (weeks)")

}
