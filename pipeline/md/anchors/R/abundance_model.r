abundance.model=function(bin.size.ifn, abundance.mat.ifn, counts.ifn, ofn.prefix)
{
  cat(sprintf("reading bin size table: %s\n", bin.size.ifn))
  sizes = read.delim(bin.size.ifn)

  cat(sprintf("reading fixed abundance matrix: %s\n", abundance.mat.ifn))
  mat = read.delim(abundance.mat.ifn)

  cat(sprintf("reading actual bin counts: %s\n", counts.ifn))
  counts = read.delim(counts.ifn)

  # total observed contacts
  N = sum(as.numeric(counts$count))

  m = merge(sizes, mat, by=c("abundance_bin1", "abundance_bin2"))

  # total potential contacts
  T = sum((as.numeric(sizes$total)))
  
  # prior contact probability (w/o abundance)
  prior = N / T

  # correction factor, so that total predicted contacts will be the observed
  factor = T / sum((as.numeric(m$total * m$probs)))

  mat$p = mat$probs * factor
  result = data.frame(abundance_bin1=mat$abundance_bin1, abundance_bin2=mat$abundance_bin2, probs=mat$p)

  ofn.f = paste(ofn.prefix, "_abundance_bin.f", sep="")
  cat(sprintf("writing result f: %s\n", ofn.f))
  write.table(result, ofn.f, quote=F, row.names=F, sep="\t")

  ofn.prior = paste(ofn.prefix, ".prior", sep="")
  cat(sprintf("writing result prior: %s\n", ofn.prior))
  write.table(prior, ofn.prior, quote=F, col.names=F, row.names=F, sep="\t")
}
