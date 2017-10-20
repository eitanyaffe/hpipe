assert=function (cond, msg="assert failed")
{
  if (!cond) {
    .Internal(stop(as.logical(T), .makeMessage(msg)) )
  }
}

# bin to a number of equal sized bins
bin.field.by.quantiles=function(values, field, nbins)
{
  cat(sprintf("binning field %s, nbins=%s\n" ,field, nbins))

  N = length(values)
  order = order(values)
  break.ind = c(1, round(1:(nbins-1)*(N/nbins)), N)

  breaks = values[order[break.ind]]
  bin.values = breaks[1:(length(breaks)-1)] + diff(breaks)/2

  bins = rep(-1, N)
  for (i in 1:nbins)
  {
    ind = order[break.ind[i]:break.ind[i+1]]
    bins[ind] = i
  }
  assert(all(bins != -1))

  return (list(bins=bins, breaks=breaks, values=bin.values, N=nbins))
}

# bin according to quantile breaks
bin.field.by.quantile.breaks=function(values, field, quantile.breaks)
{
  cat(sprintf("binning field %s, q-breaks=%s\n" ,field, paste(quantile.breaks, collapse=",")))
  nbins = length(quantile.breaks)+1

  N = length(values)
  order = order(values)
  break.ind = c(1, round(quantile.breaks*N), N)

  breaks = values[order[break.ind]]
  bin.values = breaks[1:(length(breaks)-1)] + diff(breaks)/2

  bins = rep(-1, N)
  for (i in 1:nbins)
  {
    ind = order[break.ind[i]:break.ind[i+1]]
    bins[ind] = i
  }
  assert(all(bins != -1))

  return (list(bins=bins, breaks=breaks, values=bin.values, N=nbins))
}

# bin by a vector of break values
bin.field.by.value.breaks=function(values, field, breaks)
{
  cat(sprintf("binning field %s, breaks=%s\n" ,field, paste(breaks, collapse=",")))

  bins = rep(-1, length(values))
  N = length(breaks)+1
  for (i in 1:N)
  {
    if (i > 1) from = breaks[i-1] else from = -Inf
    if (i < N) to = breaks[i] else to = Inf
    ind = (from <= values) & (values < to)
    bins[ind] = i
  }
  assert(all(bins != -1))

  breaks=c(-Inf, breaks, Inf)
  bin.values = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  bin.values[1] = min(values)
  bin.values[length(bin.values)] = max(values)

  return (list(bins=bins, breaks=breaks, values=bin.values, N=(length(breaks)-1)))
}

write.bin.features=function(prefix, table, field, bin.result, num.nas)
{
  N = bin.result$N
  bins.fn = paste(prefix, "_", field, ".bins", sep="")
  labels = character(N)
  for (i in 1:N)
  {
    from = round(bin.result$breaks[i],4)
    to = round(bin.result$breaks[i+1],4)

    if (i == 1)      labels[i] = paste("min - ", to, sep="")
    else if (i == N) labels[i] = paste(from, " - max", sep="")
    else             labels[i] = paste(from, " - ", to, sep="")
  }
  # count number of fragends in each bin
  bfield = paste(field, "bin", sep="_")
  count = integer(N)
  for (i in 1:N)
    count[i] = sum(bin.result$bins == i)

  # append NA bin
  count = c(num.nas, count)
  values = c(NA, bin.result$values)
  labels = c("NA", labels)
  result = data.frame(bin=0:N, count=count, value=values, labels=labels)
  cat(sprintf("writing table %s\n", bins.fn))
  write.table(result, bins.fn, quote=F, col.names=T, row.names=F, sep="\t")

}

bin.fields.general.inner=function(table, prefix, fields)
{
  # correlations
  if (length(fields) > 1) {
    ctable = na.omit(table[,fields])
    if (is.element("tss", fields)) ctable[,"tss"] = abs(ctable[,"tss"])
    cor = round(cor(ctable[,fields], method="spearman"),2)
    cor.fn = paste(prefix, "cor", sep=".")
    write.table(cor, cor.fn, quote=F, col.names=T, row.names=F, sep="\t")
  }

  default.nbins = 20

  for (field in fields)
  {
    values.full = table[,field]
    ind = !is.na(values.full)
    values = values.full[ind]

    #cis.breaks = c(1, 1000*1:100)
    cis.breaks = c((1:50)*2*10^3, 200000, 400000, (2:10)*10^5, (2:4)*10^6)
    pre.breaks = c((1:50)*2*10^3, 200000, 400000, (2:10)*10^5, (2:4)*10^6)

    bin.result = switch(field,
      tss_cis=bin.field.by.value.breaks(values, field, breaks=c(rev(-cis.breaks), -1, 1, cis.breaks)),
      tss_cis_active=bin.field.by.value.breaks(values, field, breaks=c(rev(-cis.breaks), -1, 1, cis.breaks)),
      random_peaks=bin.field.by.value.breaks(values, field, breaks=c(rev(-cis.breaks), cis.breaks)),
      pre=bin.field.by.value.breaks(values, field, breaks=c(rev(-pre.breaks), 0, pre.breaks)),
      tss=bin.field.by.value.breaks(values, field, breaks=c(-50000,-10000,-1,1,10000,50000)),
      bin.field.by.quantiles(values, field, nbins=default.nbins))

    # write bin results
    write.bin.features(prefix, table, field, bin.result, sum(!ind))

    # Set NA bin to zero
    bins = rep(-1, length(values.full))
    bins[which(ind)] = bin.result$bins
    bins[!ind] = 0
    assert(all(bins != -1))

    # append to table
    table = cbind(table, bins)
    colnames(table)[dim(table)[2]] = paste(field, "bin", sep="_")
  }
  return(table)
}

bin.fields.general=function(prefix="h_results/h1_cl_gc_200_ann", in.ext="fends", out.ext="binned",
                            fields=c("dnase", "K4m1", "K4m3", "ctcf"))
{
  ifn = paste(prefix, in.ext, sep=".")
  binned.fn = paste(prefix, out.ext, sep=".")

  cat(sprintf("reading table %s\n", ifn))
  table = read.delim(ifn)

  result = bin.fields.general.inner(table=table, prefix=prefix, fields=fields)

  cat(sprintf("writing table %s\n", binned.fn))
  write.table(result, binned.fn, quote=F, col.names=T, row.names=F, sep="\t")
}

###########################################################################################
# normalizes a list of fields to equal sized bins
###########################################################################################

bin.fields=function(fends.fn="h_results/h1_cl_gc_200.fends",
                    ofn.prefix="h_results/h1_cl_gc_200",
                    fields=c("frag_len", "frag_gc"), qns=c(20, 20), round=3)
{
  assert(length(fields) == length(qns))

  binned.fn = paste(ofn.prefix, ".binned", sep="")
  bin.labels.fn = paste(ofn.prefix, ".bin_labels", sep="")
  bins.fn = paste(ofn.prefix, ".bins", sep="")

  cat(sprintf("reading table %s\n", fends.fn))
  data = read.delim(fends.fn, stringsAsFactors=F)
  data = na.omit(data)

  N = dim(data)[1]
  all.labels = NULL
  mid.df = NULL
  for (i in 1:length(fields))
  {
    field = fields[i]
    qn = qns[i]

    cat(sprintf("binning field %s, qn=%s\n" ,field, qn))

    order = order(data[,field])
    breaks = c(1, round(1:(qn-1)*(N/qn)), N)

    break.values = data[order[breaks],field]
    if (!is.na(round)) break.values = round(break.values, round)

    mid.vals = break.values[1:(length(break.values)-1)] + diff(break.values)/2
    mid.df[[field]] = mid.vals

    bins = rep(-1, N)
    labels = character(qn)
    ranges = NULL
    for (i in 1:qn)
    {
      ind = order[breaks[i]:breaks[i+1]]
      bins[ind] = i

      from = break.values[i]
      to = break.values[i+1]
      ranges = rbind(ranges, c(i, from, to))
      if (i == 1)
        labels[i] = paste("min-", to, sep="")
      else if (i == qn)
        labels[i] = paste(from, "-max", sep="")
      else
        labels[i] = paste(from, "-", to, sep="")
    }
    assert(all(bins != -1))
    df = data.frame(bins)
    names(df) = paste(field, "bin", sep="_")
    data = cbind(data, df)

    df.labels = data.frame(field=paste(field, 1:qn, sep="-"), label=labels)

    all.labels = rbind(all.labels, df.labels)

    # write ranges
    ranges.ofn = paste(ofn.prefix, "_",field, ".bin_ranges", sep="")
    colnames(ranges) = c("bin", "start", "end")
    write.table(ranges, ranges.ofn, quote=F, col.names=T, row.names=F, sep="\t")
  }

  cat(sprintf("writing table %s\n", binned.fn))
  write.table(data, binned.fn, quote=F, col.names=T, row.names=F, sep="\t")

  # write labels
  cat(sprintf("writing table %s\n", bin.labels.fn))
  write.table(all.labels, bin.labels.fn, quote=F, col.names=T, row.names=F, sep="\t")

  # count number of fragend in each bin
  bfields = paste(fields, "bin", sep="_")
  keys = unique(data[,bfields])
  qn.ranges = list()
  for (i in 1:length(fields))
    qn.ranges[[i]] = 1:qns[i]
  keys = expand.grid(qn.ranges)
  M = dim(keys)[1]

  count = rep(round(N/M), M)

  bins = cbind(keys, count)
  names(bins) = c(bfields, "count")
  bins = cbind(bins, as.data.frame(mid.df))

  cat(sprintf("writing table %s\n", bins.fn))
  write.table(bins, bins.fn, quote=F, col.names=T, row.names=F, sep="\t")
}

bin.fields.mfn=function(
    fends.fn="h_results/h1_cl_gc_200.fends",
    ofn.prefix="h_results/h1_cl_gc_200",
    mfn, round=3)
{
    mtable = read.delim(mfn, stringsAsFactors=F)
    fields = mtable[, "field"]

    ofn = paste(ofn.prefix, "binned", sep=".")
    if (any(mtable$type == "optimize") > 0) {
        rfields = mtable[mtable$type == "optimize", "raw_field"]
        qns = mtable[mtable$type == "optimize", "size"]
        ix = qns>0
        bin.fields(fends.fn=fends.fn, ofn.prefix=ofn.prefix, fields=rfields[ix], qns=qns[ix], round=round)
    } else {
        cat(sprintf("creating binned table: %s\n", ofn))
        system(paste("cp -r", fends.fn, ofn))
    }
}
