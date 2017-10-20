plot.gs=function(ifn.ref.prefix, ifn.anchor.prefix, min.k, max.k, aai.breaks, fdir)
{
    for (k in min.k:max.k) {
        ifn.ref = paste(ifn.ref.prefix, ".", k, sep="")
        ifn.anchor = paste(ifn.anchor.prefix, ".", k, sep="")
        df.ref = load.table(ifn.ref)
        df.anchor = load.table(ifn.anchor)
    }
}
