library(data.table)
library(stringi)

wd <- getwd()
p <- fread(file.path(wd, "matFetFrags", "informativeReads.pileup"))
setnames(p, c("chr", "pos", "ref", "dep", "pileup", "qual", "reads"))
b <- fread(file.path(wd, "matFetFrags", "informativeSites.bed"))
b[ , V2 := NULL]
setnames(b, c("chr", "pos", "keep"))
setkey(p, chr, pos)
setkey(b, chr, pos)
p <- b[p]
## Remove loci that had indels for simplicity
p <- p[!stri_detect_regex(pileup, "\\+|\\-|\\*")]

## Figure out how to correctly identify read positions
# p[ , np := nchar(pileup)]
# p[ , nd := stri_count_regex(pileup, "\\$")]
# p[ , nc := stri_count_regex(pileup, "\\^")]
# p[ , nr := sapply(stri_split_regex(reads, "\\,"), length)]
# p[ , all(np == nd + nc*2 + nr)]

p[ , simplePileup := stri_replace_all_regex(pileup, "\\$|\\^.", "")]
# p[ , all(nchar(simplePileup) == nr)]

p[ , rgx := ifelse(grepl("Ref", keep), "\\.|\\,", "[AaCcTtGg]")]
p[ , sel := stri_locate_all_regex(simplePileup, rgx)]
p[ , sel := lapply(sel, function(x) x[ , "start"])]
p[ , sub := mapply("[", stri_split_regex(reads, "\\,"), sel)]
p[ , smp := ifelse(grepl("fet", keep), "fet", "mat")]

fetRds <- p[smp == "fet", unique(unlist(sub))]
matRds <- p[smp == "mat", unique(unlist(sub))]

## Remove the overlapping reads
length(intersect(fetRds, matRds))
fetRds <- setdiff(fetRds, matRds)
matRds <- setdiff(matRds, fetRds)

## Write files to pass to picard tools
wrtRds <- function(x, nm) {
  write.table(x, nm, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
wrtRds(fetRds, "fetalReads.txt")
wrtRds(matRds, "maternalReads.txt")



