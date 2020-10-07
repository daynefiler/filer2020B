#!/usr/bin/env Rscript

"Usage: getAlleleCounts.R -v <vcf> -o <outfile> 

Options:
  -h --help
  -v <vcf> --vcf <vcf> VCF to extract allele depths from
  -o <outfile> --outfile <outfile>  output filepath
  
Notes: 
  Output saves as RDS file with a data.table object.
" -> doc

suppressPackageStartupMessages({
  stopifnot(require(docopt, quietly = TRUE))
  opt <- docopt::docopt(doc)
  stopifnot(require(data.table, quietly = TRUE))
  stopifnot(require(VariantAnnotation, quietly = TRUE))
})

file.create(opt[["--outfile"]])

d <- readVcf(opt[["--vcf"]])

if (ncol(geno(d)$AD) > 1) stop("vcf should only contain 1 sample")
if (max(lengths(geno(d)$AD)) > 2) {
  stop("vcf should be normed with multi-alleleic sites ",
       "split into multiple records")
}

a <- as.data.table(do.call(rbind, geno(d)$AD[ , 1]), keep.rownames = TRUE)
setnames(a, c("id", "ref", "alt"))
a[ , gt := as.vector(geno(d)$GT)]
a[ , paste0("pl", 1:3) := as.data.table(do.call(rbind, geno(d)$PL[ , 1]))]

saveRDS(a, file = opt[["--outfile"]])






