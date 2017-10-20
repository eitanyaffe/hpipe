coverage.table=function(table, idir, ofn)
{
  cat(sprintf("reading contig table: %s\n", table))
  t = read.delim(table)

  cat(sprintf("going over %d contig coverage files\n", dim(t)[1]))
  for (i in 1:dim(t)[1]) {
    id = t[i,"contig"]
    ifn = paste(idir, "/", id, sep="")
    x = read.delim(ifn)
    t[i,"sum"] = sum(x$count)
    t[i,"median"] = median(x$count)
    t[i,"sd"] = sd(x$count)
  }
  t$reads.per.bp = t$sum / t$length
  t$abundance = t$reads.per.bp / sum(t$reads.per.bp)
  t$relative.abundance = t$abundance / median(t$abundance)

  total.length = sum(t$length)
  total.reads = sum(t$sum)
  t$abundance.enrichment = log10((1 + t$sum) / (1 + (total.reads * (t$length / total.length))))

  fields = c("contig", "length", "median", "sd", "sum", "reads.per.bp", "abundance", "relative.abundance",  "abundance.enrichment")
  cat(sprintf("writing contig coverage summary table: %s\n", ofn))
  write.table(t[,fields], ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

cell.coverage.table=function(table, idir, ofn.median, ofn.sum)
{
  options(stringsAsFactors=F)
  cat(sprintf("reading contig table: %s\n", table))

  result = read.delim(table)
  result = result[,c("contig", "length")]
  result = result[1:1000,]

  rmed = result
  rsum = result

  cat(sprintf("going over %d contigs ...\n", dim(result)[1]))
  for (i in 1:dim(result)[1]) {
    id = result[i,"contig"]
    ifn = paste(idir, "/", id, sep="")
    x = read.delim(ifn)
    n = dim(x)[2] - 1
    for (j in 1:n) {
      field = paste("cell", j, sep="")
      rmed[i,field] = median(x[,field])
      rsum[i,field] = sum(x[,field])
    }
  }
  cat(sprintf("writing contig coverage summary median table: %s\n", ofn.median))
  write.table(rmed, ofn.median, quote=F, col.names=T, row.names=F, sep="\t")

  cat(sprintf("writing contig coverage summary sum table: %s\n", ofn.sum))
  write.table(rsum, ofn.sum, quote=F, col.names=T, row.names=F, sep="\t")

}
