library(data.table)
wd <- getwd() ## Snakemake directory

sampleMeta <- data.table(smp = c("S1", "S2", sprintf("FES-0034-%d", c(0:2, 4))))
sampleMeta[ , case := c(1, 2, rep(3, 4))]
sampleMeta[ ,
            source := c("cell-free", "cell-free", "newborn", "maternal",
                        "paternal", "cell-free")]
save(sampleMeta, file = file.path(wd, "forLetter", "sampleMeta.RData"))

##----------------------------------------------------------------------------##
## readSmry
##----------------------------------------------------------------------------##

## Read in markdup files
mdFls <- list.files(file.path(wd, "markdup"), full.names = TRUE)
readMarkDup <- function(f) {
  d <- as.data.table(read.table(f, nrows = 1, header = TRUE))
  d[ , f := basename(f)]
  d[]
}
md <- rbindlist(lapply(mdFls, readMarkDup))
md[ , smp := sub(".markdup.+$", "", f)]

## Read in align metrics files
amFls <- Sys.glob(file.path(wd, "*", "*.alignMetrics"))
readAlign <- function(f) {
  d <- fread(f)
  d[ , f := basename(f)]
  d[]
}
am <- rbindlist(lapply(amFls, readAlign))
am[ , smp := sub(".keep.+$|.recal.+$", "", f)]
am[ , set := ifelse(grepl("recal", f), "recal", "keep")]

## Create summary object
rp <- md[ , .(smp, ttl = READ_PAIRS_EXAMINED*2, pctDup = PERCENT_DUPLICATION)]
rp[ , dup := ttl*pctDup]
setkey(rp, smp)
setkey(am, smp)
rp <- am[set == "keep" & CATEGORY == "PAIR", .(smp, keep = TOTAL_READS)][rp]
rp[ , pctFlt := (ttl - dup - keep)/ttl*100]
rp[ , pctDup := pctDup *100]

readSmry <- list(summary = rp,
                 markDuplicates = md,
                 alignMetrics = am)

save(readSmry, file = file.path(wd, "forLetter", "readSmry.RData"))

##----------------------------------------------------------------------------##
## gt
##----------------------------------------------------------------------------##

gtFls <- list.files(file.path(wd, "final"), full.names = TRUE)
rdFinal <- function(fl) {
  x <- readRDS(fl)
  x[ , smp := sub(".final", "", basename(fl))]
}

gt <- rbindlist(lapply(gtFls, rdFinal))
gt <- gt[!is.na(adep)]
gt[ , c("chr", "var") := tstrsplit(id, ":")]
gt[ , chr := as.integer(gsub("^NC_0+|[:.:][0-9][0-9]$", "", chr))]
gt[ , c("pos", "bp") := tstrsplit(var, "_", type.convert = TRUE)]
gt[ , pos := round(pos, -4)]
gt[ , c("var", "bp") := NULL]
gt[ , varid := .GRP, by = id]
rs140468248 <- gt[id == "NC_000001.11:42755598_C/A", varid]
gt[ , id := NULL]
setkey(gt, smp, chr, pos, varid)
setcolorder(gt)

save(gt, file = file.path(wd, "forLetter", "gt.RData"))
save(rs140468248, file = file.path(wd, "forLetter", "rs140468248.RData"))
