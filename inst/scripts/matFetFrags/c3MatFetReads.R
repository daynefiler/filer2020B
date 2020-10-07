library(Rsamtools)
library(data.table)

rdBam <- function(fl) {
  p <- ScanBamParam(what = c("isize", "qname", "rname", "pos", "mapq"))
  x <- as.data.table(scanBam(fl, param = p)[[1]])
  x <- x[ ,
         .(pos = min(pos), mapq = mean(mapq), isize = max(isize)),
         by = .(qname, rname)]
  x[]
}

wd <- getwd()
mat <- rdBam(file.path(wd, "matFetFrags", "maternalReads.bam"))
mat[ , source := "maternal"]
fet <- rdBam(file.path(wd, "matFetFrags", "fetalReads.bam"))
fet[ , source := "fetal"]
c3MatFetReads <- rbind(mat, fet)
save(c3MatFetReads, file = file.path(wd, "forLetter", "c3MatFetReads.RData"))
