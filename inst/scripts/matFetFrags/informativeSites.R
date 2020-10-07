library(data.table)
wd <- getwd() ## Snakemake directory

gtFls <- file.path(wd, "final", sprintf("FES-0034-%d.final", 0:1))
rdFinal <- function(fl) {
  x <- readRDS(fl)
  x[ , smp := sub(".final", "", basename(fl))]
}

gt <- rbindlist(lapply(gtFls, rdFinal))
gt <- gt[!is.na(adep)]
gt[ , c("chr", "var") := tstrsplit(id, ":")]
gt[ , c("pos", "bp") := tstrsplit(var, "_", type.convert = TRUE)]

setorder(gt, smp, chr, pos, -alt)
myfirst <- function(x) x[1]
tbl <- dcast(gt, chr + pos ~ smp, value.var = "gt", fun.aggregate = myfirst)
setnames(tbl, c("chr", "pos", "fet", "mat"))
tbl[ , keep := NA_character_]
tbl[fet == "0/1" & is.na(mat),   keep := "fetAlt"]
tbl[fet == "0/1" & mat == "0/0", keep := "fetAlt"]
tbl[fet == "0/1" & mat == "1/1", keep := "fetRef"]
tbl[mat == "0/1" & is.na(fet),   keep := "matAlt"]
tbl[mat == "0/1" & fet == "0/0", keep := "matAlt"]
tbl[mat == "0/1" & fet == "1/1", keep := "matRef"]

md <- fread(file.path(wd, "depth", "FES-0034-1.all.depth"))
fd <- fread(file.path(wd, "depth", "FES-0034-0.all.depth"))
setnames(md, c("chr", "pos", "matDep"))
setnames(fd, c("chr", "pos", "fetDep"))
setkey(md, chr, pos)
setkey(fd, chr, pos)

tbl <- tbl[!is.na(keep)]
tbl <- md[tbl]
tbl <- fd[tbl]
tbl <- tbl[fetDep > 20 & matDep > 20]

write.table(tbl[ , .(chr, pos - 1, pos, keep)],
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            file = file.path(wd, "matFetFrags", "informativeSites.bed"))

