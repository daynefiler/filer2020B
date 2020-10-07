#!/usr/bin/env Rscript

"Usage: getAdepFdist.R -a <adep>  -f <fdist>  [-t <tdep>] [-m <mdep>] -o <outfile>

Options:
  -h --help
  -a <adep> --adep <adep> allele dep file 
  -f <fdist> --fdist <fdist> fragment length distribution file
  -t <tdep> --tdep <tdep> the minimum total depth for inclusion
  -m <mdep> --mdep <mdep> the minimum depth for both alt & ref for inclusion
  -o <outfile> --outfile <outfile> the output file
  
Notes: 
  Output saves as RDS file with a data.table object.
" -> doc

suppressPackageStartupMessages({
  stopifnot(require(docopt, quietly = TRUE))
  opt <- docopt::docopt(doc)
  stopifnot(require(data.table, quietly = TRUE))
})

file.create(opt[["--outfile"]])

f <- readRDS(opt[["--fdist"]])
a <- readRDS(opt[["--adep"]])
a <- a[ref + alt >= opt[["--tdep"]] & 
         ref >= opt[["--mdep"]] & 
           alt >= opt[["--mdep"]]]
a[ , c("chr", "pos") := tstrsplit(id, split = ":")]
a[ , c("pos") := tstrsplit(pos, split = "_", type.convert = TRUE, keep = 1)]

setkey(a, chr, pos)
setkey(f, chr, pos)
a <- f[a]
set(a, j = "chr", value = NULL)
set(a, j = "pos", value = NULL)
setkey(a, id)
setcolorder(a)

saveRDS(a, file = opt[["--outfile"]])


