bin.abundance=function(fends.ifn, fends.ofn, bin.ranges.ofn, correction.ofn, mfn)
{
  cat(sprintf("reading fends table: %s\n", fends.ifn))
  table = read.delim(fends.ifn)


  mfile = read.delim(mfn)
  im = match("abundance", mfile$raw_field)
  if (is.na(im))
      stop("abundance field not found in model file")
  n.bins = mfile[im, "size"]
  if (!is.numeric(n.bins))
      stop("n bins not defined")

  cat(sprintf("number of bins: %d\n", n.bins))
  N = n.bins

  rr = range(table$abundance)
  breaks = seq(rr[1], rr[2], length.out=N+1)
  bins = data.frame(bin=1:N, start=breaks[1:N], end=breaks[2:(N+1)])
  middle.log = (bins$start + bins$end) / 2
  middle = 10^middle.log

  cmat = expand.grid(1:N, 1:N)
  names(cmat) = c("abundance_bin1", "abundance_bin2")
  factors = expand.grid(middle, middle)
  cmat$probs = factors[,1] * factors[,2]
  cmat = cmat[cmat[,1] <= cmat[,2],]

  # make factors distribute around 1
  cmat$probs = cmat$probs / median(cmat$probs)

  cat(sprintf("saving correction matrix file: %s\n", correction.ofn))
  write.table(cmat, correction.ofn, quote=F, col.names=T, row.names=F, sep="\t")

  table$abundance_bin = as.numeric(cut(table$abundance, breaks, include.lowest=T))
  if (any(is.na(table$abundance_bin)))
    stop("NA found in abundance_bin")
  cat(sprintf("saving fends file: %s\n", fends.ofn))
  write.table(table, fends.ofn, quote=F, col.names=T, row.names=F, sep="\t")

  cat(sprintf("saving bin range file: %s\n", bin.ranges.ofn))
  bins$start = round(10^bins$start, 4)
  bins$end = round(10^bins$end, 4)
  write.table(bins, bin.ranges.ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
