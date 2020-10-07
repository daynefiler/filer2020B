#!/usr/bin/env Rscript

"Usage: getFragDist.R -a <adep> -l <ldep> -s <sdep> [-r <region>] -o <outfile> 

Options:
  -h --help
  -a <adep> --alldep <adep> depth file for all fragments
  -l <ldep> --longdep <ldep> depth file for the long fragments
  -s <sdep> --shortdep <sdep> depth file for the short fragments
  -r <region> --region <region> BED file to calculate depths over
  -o <outfile> --outfile <outfile>  output filepath
  
Notes: 
  Output saves as RDS file with a data.table object. Omitting --region
  will default to using the intersect of all three depth files. For file
  latency purposes, the program starts by creating an empty file.
  Requires R, 'data.table' (R package), and 'docopt' (R package) to be 
  installed Depth files are expected to be 0-based coordinates. BED files
  for the --region parameter are expected to be 0-based for the start and
  1-based for the end The output is given in 0-based coordinates.
" -> doc

suppressPackageStartupMessages({
  stopifnot(require(docopt, quietly = TRUE))
  opt <- docopt::docopt(doc)
  stopifnot(require(data.table, quietly = TRUE))
})

file.create(opt[["--outfile"]])

cnms <- c("chr", "pos")

adep <- fread(opt[["--alldep"]],   col.names = c(cnms, "adep"), key = cnms)
sdep <- fread(opt[["--shortdep"]], col.names = c(cnms, "sdep"), key = cnms)
ldep <- fread(opt[["--longdep"]],  col.names = c(cnms, "ldep"), key = cnms)

if (is.null(opt[["--region"]])) {
  dep <- merge(ldep, merge(adep, sdep, all = TRUE), all = TRUE)
} else {
  dep <- fread(opt[["--region"]], 
               col.names = c("chr", "start", "end"), 
               select = 1:3)
  dep <- dep[ , .(chr = chr, pos = seq(start, end - 1)), by = 1:nrow(dep)]
  dep[ , nrow := NULL]
  setkeyv(dep, cols = cnms)
  dep <- merge(dep, adep, all.x = TRUE)
  dep <- merge(dep, sdep, all.x = TRUE)
  dep <- merge(dep, ldep, all.x = TRUE)
}

dep[is.na(adep), adep := 0]
dep[is.na(sdep), sdep := 0]
dep[is.na(ldep), ldep := 0]

saveRDS(dep, file = opt[["--outfile"]])






